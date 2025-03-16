from collections import defaultdict
import sys
#chromosome position strand sequence_context mc cov methylated
#chr1	3001277	+	CGA	65	66	1
site_dict = defaultdict(list)
with open(sys.argv[1]) as fin:
    for line in fin:
        linL = line.strip('\n').split('\t')
        chrom = linL[0]
        site = int(linL[1])
        strand = linL[2]
        mc = int(linL[4])
        cov = int(linL[5])
        if strand=='-':
            m_site = chrom + '_' + str(site-1)
            if m_site in site_dict:
                site_dict[m_site][0] += mc
                site_dict[m_site][1] += cov
            else:
                site_dict[m_site].append(mc)
                site_dict[m_site].append(cov)
        else:
            m_site = chrom + '_' + str(site)
            if m_site in site_dict:
                site_dict[m_site][0] += mc
                site_dict[m_site][1] += cov
            else:
                site_dict[m_site].append(mc)
                site_dict[m_site].append(cov)
fo = open(sys.argv[2], 'w')
for ele in site_dict:
    chrom, start = ele.split('_')[0],  ele.split('_')[1]
    mc = site_dict[ele][0]
    cov = site_dict[ele][1]
    fo.write(chrom + '\t' + str(start) + '\t' + str(start) + '\t' + str(mc) + '\t' + str(cov))
    fo.write('\n')
