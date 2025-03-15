from collections import defaultdict
import re
import pandas as pd
import pysam
import sys
import argparse
import multiprocess
import time
import subprocess
import signal
import os
import traceback


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


args = sys.argv[1:]
parser = argparse.ArgumentParser(
    description="The code to filter reads for RNA after STAR align\
        (from snmCAT-seq(mct_pip)),\
    Filtered RNA reads containing full mCH, \
    speed up by processing chromosome in parallel (according to species), \
    Update in 2023.09.12 by Yingchuo Hu. \
    This code is further update by Liuhao Ren at 240809 \
    for restarting apply_async at Exception.\
    Last editting time: 2024-08-13 14:14:47 .",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i", "--sampleid",
                    required=True,
                    help="The Sampleid.")

parser.add_argument("-s", "--species",
                    required=True,
                    help="The species(hg38/mm10/lambda/puc19/clai).")

parser.add_argument("-b", "--input_bam",
                    required=True,
                    help="The input bam (STAR align).")

parser.add_argument("-o", "--output_bam",
                    type=str,
                    required=True,
                    help="The output bam(filtered RNA reads).")

parser.add_argument("-m", "--mC_min",
                    type = float,
                    default=0.9,
                    help="mc_rate_min_threshold")

parser.add_argument("-c", "--cov_min",
                    type = int,
                    default=3,
                    help="cov_min_threshold")

parser.add_argument("--sequence",
                    type = str,
                    default="PE",
                    help="SE/PE")

parser.add_argument("-t", "--thread",
                    type = int,
                    default=10,
                    help="thread for split chromosome")


args = parser.parse_args()
sampleid = args.sampleid
species = args.species
input_bam = args.input_bam
output_bam = args.output_bam
mc_rate_min_threshold = args.mC_min
cov_min_threshold = args.cov_min
sequence = args.sequence
threads = args.thread

output_txt = output_bam + ".reads_CH.txt"

REVERSE_READ_MCH_CONTEXT = {'CA', 'CC', 'CT'}
FORWARD_READ_MCH_CONTEXT = {'AG', 'TG', 'GG'}

spikeIN={"lambda","fullpuc19","clai"}

'''
For RNA library (STAR align)

species: 
hg38,mm10: by chrom (Full mCH/hmCH)

for 5hmC :
lambda,puc19,clai: (Full hmCH)

for 5mC:
lambda,puc19,clai: (Full mCH)

'''


def single_read_mch_level(read):

    # ref seq is parsed based on read seq and MD tag, and do not depend on reverse or not
    
    ref_seq = read.get_reference_sequence().upper()
    ref_pos = read.get_reference_positions()
    # use dict instead of string is because ref_seq could contain blocks when skip happen
    ref_seq_dict = {pos: base for pos, base in zip(ref_pos, ref_seq)}
    read_seq = read.seq.upper()

    # only count mCH
    mch = 0
    cov = 0
    other_snp = 0
    if read.is_reverse:  # read in reverse strand
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ##ref_read_pair = ref_base + read_base
            ref_read_pair = ref_base.upper() + read_base
            try:
                ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos + 1]
                if ref_context not in REVERSE_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue
            if ref_read_pair == 'CC':  # C to C means unconverted and methylated
                cov += 1
                mch += 1
            elif ref_read_pair == 'CT':  # C to T means converted and un-methylated
                cov += 1
            else:
                # other kinds of SNPs, do not count to cov
                other_snp += 1
                pass
    else:  # read in forward strand
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos - 1] + ref_seq_dict[ref_pos]
                if ref_context not in FORWARD_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue
            if ref_read_pair == 'GG':  # G to G means unconverted and methylated
                cov += 1
                mch += 1
            elif ref_read_pair == 'GA':  # G to A means converted and un-methylated
                cov += 1
            else:
                # other kinds of SNPs, do not count to cov
                other_snp += 1
                pass

    read_mch_frac = (mch / cov) if cov > 0 else 0
    ##print(read.query_name, str(mch), str(mch),str(cov), str(read_mch_frac))
    return read_mch_frac, mch, cov, other_snp 



def select_rna_reads_normal(input_bam,
                           chrom,
                           bam_out_tmp,
                           txt_out_tmp):
    
    global lock
    lock = multiprocess.Lock()
    

    reads = 0
    Filtered_reads = 0
    with pysam.AlignmentFile(input_bam) as bam:
        with pysam.AlignmentFile(bam_out_tmp, header=bam.header,
                                 mode='wb') as out_bam, \
        open(txt_out_tmp, "w") as outfile:
            for read in bam.fetch(chrom):
                reads += 1
                read_mch_rate, mch, cov, other_snp = single_read_mch_level(read)
                read_id = read.query_name

                ##print(read.query_name, str(read_mch_rate),str(cov))
                
                lock.acquire()
                outfile.write(f'{read_id}\t{read_mch_rate}\t{mch}\t{cov}\n')
                lock.release()
                
                # split reads
                if (read_mch_rate <
                    mc_rate_min_threshold) or (cov < cov_min_threshold):
                    continue
                    
                lock.acquire()
                Filtered_reads += 1
                out_bam.write(read)
                lock.release()
                 
        
    return


def read_methyl (read, type) :
    ref_seq = read.get_reference_sequence().upper()
    ref_pos = read.get_reference_positions()
    ref_seq_dict = {pos: base for pos, base in zip(ref_pos, ref_seq)}
    read_seq = read.seq.upper()
    read_id = read.query_name # assigned but never used
    CH_ref_pos = defaultdict()

    # only count mCH

    if type == "CtoT":  # read CtoT conversion

        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ##ref_read_pair = ref_base + read_base
            ref_read_pair = ref_base.upper() + read_base
            try:
                ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos + 1]
                
                if ref_context not in REVERSE_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue

            if ref_read_pair == 'CC':  # C to C means unconverted and methylated
                CH_ref_pos[ref_pos]= "m" 
            elif ref_read_pair == 'CT':  # C to T means converted and un-methylated
                CH_ref_pos[ref_pos]= "u" 
            else:
                # other kinds of SNPs, do not count to cov
                CH_ref_pos[ref_pos] = read_base
                pass
            
        return(CH_ref_pos)

    else:
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos - 1] + ref_seq_dict[ref_pos]
                if ref_context not in FORWARD_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue
            if ref_read_pair == 'GG':  # G to G means unconverted and methylated
                CH_ref_pos[ref_pos]= "m" 
            elif ref_read_pair == 'GA':  # G to A means converted and un-methylated
                CH_ref_pos[ref_pos]= "u" 
            else:
                # other kinds of SNPs, do not count to cov
                CH_ref_pos[ref_pos]= read_base
                pass

        return(CH_ref_pos)

def select_rna_reads_PE(input_bam,
                           chrom,
                           bam_out_tmp,
                           txt_out_tmp):
    
    global lock
    lock = multiprocess.Lock()
    
    read1 = 0
    read2 = 0
    read_mate_unmap = 0
    read1_proper_pairs = 0
    Filtered_reads = 0
    with pysam.AlignmentFile(input_bam) as bam:
        with pysam.AlignmentFile(bam_out_tmp, header=bam.header,
                                 mode='wb') as out_bam, \
        open(txt_out_tmp, "w") as outfile:
            for read in bam.fetch(chrom):
                if read.is_read1:
                    read1 +=1
                    if read.mate_is_unmapped:
                        try:
                            read_mate_unmap+=1
                            read_id=read.query_name
                            outfile.write(f'{read_id},"mate_unmap"\n')
                        except:
                            pass                        
                        continue
                    try:
                        read1_proper_pairs += 1
                        read_id = read.query_name
                        mate = bam.mate(read)
                        ##print(read.query_name)
                        if read.is_reverse :
                            read1_CH = read_methyl(read,"CtoT")
                            read2_CH = read_methyl(mate,"GtoA")
                            pair_CH = read1_CH
                        else :
                            read1_CH = read_methyl(read,"GtoA")
                            read2_CH = read_methyl(mate,"CtoT")
                            pair_CH = read1_CH
                        for ref_pos, state in read2_CH.items():
                            if ref_pos in read1_CH :
                                if state == read1_CH[ref_pos] and state == "m" :
                                    continue
                                elif state == read1_CH[ref_pos] and state == "u" :
                                    continue
                                else:
                                    del pair_CH[ref_pos]
                            else:
                                pair_CH[ref_pos] =state

                        mch = "".join(list(pair_CH.values())).count("m")
                        uch = "".join(list(pair_CH.values())).count("u") 
                        cov = mch + uch
                        read_mch_rate = (mch / cov) if cov > 0 else 0
                        
                        lock.acquire()
                        outfile.write(f'{read_id}\t{read_mch_rate}\t{mch}\t{cov}\n')
                        lock.release()
                        
                        if (read_mch_rate < mc_rate_min_threshold) or (cov < cov_min_threshold):
                            continue
                        Filtered_reads += 2
                        
                        lock.acquire()
                        out_bam.write(read)
                        out_bam.write(mate)
                        lock.release()
                        
                    except ValueError as e:
                        print(read.query_name)
                        print(f"Error: {e}")
                        try:
                            lock.acquire()
                            outfile.write(f'{read_id}\t"mate_unfind"\n')
                            lock.release()
                        except:
                            pass
                        continue

                else:
                    read2 +=1
                    if read.mate_is_unmapped:
                        try:
                            read_mate_unmap+=1
                            read_id=read.query_name
                            lock.acquire()
                            outfile.write(f'{read_id}\t"mate_unmap"\n')
                            lock.release()
                        except:
                            pass
                        continue
                    try:
                        mate = bam.mate(read)
                    except ValueError as e:
                        print(read.query_name)
                        print(f"Error: {e}")
                        try:
                            read_id=read.query_name                        
                            lock.acquire()
                            outfile.write(f'{read_id}\t"mate_unfind"\n')
                            lock.release()
                        except:
                            pass
                        continue
                    
            
            
    print(f"chromosome {chrom} done successfully, with {read1} read1, {read2} read2, {read_mate_unmap} read_mate_unmap, {read1_proper_pairs} read1_proper_paris.\n")
    return


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable=callable
        return
    def __call__(self,*args,**kwargs):
        try:
            result=self.__callable(*args,**kwargs)
        except Exception:
            print(traceback.format_exc())
        return result
    pass


def run_all():

    ##t0=clock()
    t0=time.perf_counter()
    
    samfile = pysam.AlignmentFile(input_bam, "rb")
    assert samfile.check_index(),  "ErrorType: %s file does not have index file." % input_bam
    print ('Start Analysis threads: %s' % threads)

    print("Input is " + input_bam,
              "Ouput is " + output_bam,
              "Species is " + species,
              "mc_rate_min_threshold is " + str(mc_rate_min_threshold), 
              "cov_min_threshold is " + str(cov_min_threshold),
              "Library type is " + sequence)
            
    
    header = pd.DataFrame(samfile.header['SQ'])
    Chrs = header['SN'].tolist()
    samfile.close()

    #multiprocess.log_to_stderr() #Some delay in stderr causing reruning of successful subprocess?

    while Chrs:
        pool = multiprocess.Pool(threads,init_worker)
        result_dict = {}
        tmp_ReadOut = set()
        tmp_Readtxt = set()
        for chrom in Chrs[:]:
            ###print("all chrom is " + chrom)
            if species not in spikeIN:
                if not re.search(r'^chr',chrom):
                    Chrs.remove(chrom)
                    continue        
            if chrom == "chrM":
                Chrs.remove(chrom)
                continue
                
            bam_out_tmp = output_bam + '.' + chrom + ".bam"
            txt_out_tmp = output_txt + '.' + chrom + ".txt"
            tmp_ReadOut.add(bam_out_tmp)
            tmp_Readtxt.add(txt_out_tmp)
            
            # Only select_rna_reads_PE is fully implemented.
            if sequence == 'SE':
                result_dict[chrom] = pool.apply_async(func=LogExceptions(select_rna_reads_SE), args=(input_bam, chrom, bam_out_tmp, txt_out_tmp))
            elif sequence == 'PE':
                result_dict[chrom] = pool.apply_async(func=LogExceptions(select_rna_reads_PE), args=(input_bam, chrom, bam_out_tmp, txt_out_tmp))
        
        pool.close()
        pool.join()
        for chrom in Chrs[:]:
            if result_dict[chrom].successful():
                Chrs.remove(chrom)
        print(f"{Chrs} remaining to be dealt.\n")    
    
        
        
    try:
        print ("Waiting 1 seconds")
        time.sleep(1)
    except KeyboardInterrupt:
        print ("Caught KeyboardInterrupt, terminating workers")
    else:
        print ("Quitting normally")
        print ('Pool finished')

    
    ###print("tmp_ReadOut: ")
    ###print(tmp_ReadOut)

    ##t1=clock()
    start_time = time.perf_counter()
    print ("Filter Bam file by mCHs spend %s" % str(start_time-t0))

    Args_m = ['samtools merge -f -@ %s %s %s' % (threads, output_bam, " ".join(tmp_ReadOut))]
    Args_txt = ['cat %s > %s' % (" ".join(tmp_Readtxt) , output_txt)]
    subprocess.check_call(Args_m, shell=True)
    subprocess.check_call(Args_txt, shell=True)
         
    for temp_file in tmp_ReadOut:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    for temp_file in tmp_Readtxt:
        if os.path.exists(temp_file):
            os.remove(temp_file)      
    
    ##t2=clock()
    end_time = time.perf_counter()
    
    elapsed_time = end_time - start_time

    print ("Final merged Bam file %s in %s" % (output_bam, str(elapsed_time)))


if __name__ == "__main__":

    run_all()


