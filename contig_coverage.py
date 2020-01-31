#!/usr/bin/env python3

import sys
from pathlib import Path


def get_coverage(fname):
    min_contig_size = 1000
    contig_count = 0
    total_coverage = 0.
    with fname.open("r") as f:
        next(f)
        for line in f:
            cols = line.split('\t')
            contig_size = int(cols[0].split('_')[3])
            if contig_size < min_contig_size:
                continue
            contig_count += 1
            total_coverage += float(cols[2])
    print(f"{total_coverage/contig_count*100:.2f}%")


if __name__ == "__main__":
    get_coverage(Path(sys.argv[1]))
