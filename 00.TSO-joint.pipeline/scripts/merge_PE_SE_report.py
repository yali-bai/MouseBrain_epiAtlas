# -*- coding: utf-8 -*-


import sys
import os
import time
if len(sys.argv) < 4:
    print(sys.argv[0] + ' PE_report SE_report merge_report\n')
    exit()
PE_report = sys.argv[1]
SE_report = sys.argv[2]
start = time.time()
merge_report = open(sys.argv[3],'w')
row_num = os.popen('wc -l {}'.format(PE_report)).read()
row_num = int(row_num.split()[0])
print('report file row num:',row_num)
pe_file = open(PE_report,'r')
se_file = open(SE_report,'r')
for i in range(row_num):
    line_pe = pe_file.readline().strip().split('\t')
    line_se = se_file.readline().strip().split('\t')
    chr1,pos1,strand1,methyl1,unmethyl1,CG1,CGN1 = line_pe[:]
    chr2,pos2,strand2,methyl2,unmethyl2,CG2,CGN2 = line_se[:]
    key1 = chr1 + '_' + pos1 + '_' + strand1 + '_' + CG1 + '_' + CGN1
    key2 = chr2 + '_' + pos2 + '_' + strand2 + '_' + CG2 + '_' + CGN2
    assert key1 == key2
    methy = int(methyl1) + int(methyl2)
    unmethyl = int(unmethyl1)+ int(unmethyl2)
    outline = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr1, pos1, strand1, methy, unmethyl, CG1, CGN1)
    merge_report.write(outline)
pe_file.close()
se_file.close()
merge_report.close()
end = time.time()
print('merge use time:',end-start)

