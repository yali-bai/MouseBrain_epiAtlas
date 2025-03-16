import pyBigWig
import os
import sys

def process_bw(input_path, output_path):
    bw_in = pyBigWig.open(input_path)
    chroms = bw_in.chroms()
    bw_out = pyBigWig.open(output_path, "w")
    bw_out.addHeader(list(chroms.items()))

    for chrom in chroms:
        intervals = bw_in.intervals(chrom)
        if intervals is None:
            continue
        all_starts = []
        all_values = []
        for start, end, value in intervals:
            starts = list(range(start, end))
            values = [value] * (end - start)
            all_starts.extend(starts)
            all_values.extend(values)
        if all_starts:
            bw_out.addEntries(chrom, all_starts, values=all_values, span=1)

    bw_in.close()
    bw_out.close()
    print(f"Processed {input_path} and saved to {output_path}")

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_bw(input_file, output_file)




