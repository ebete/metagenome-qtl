#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


def get_coverage_per_gene(bed_file):
    print("orf", "avg_depth", "coverage", sep="\t")

    with bed_file.open("r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")

        covered = 0
        total = 0
        gene_nucl = 0
        last_gene = ""
        for record in handle:
            cur_gene = record[0]
            start_idx = int(record[1])
            end_idx = int(record[2])
            span = end_idx - start_idx
            read_depth = int(record[3])

            if cur_gene != last_gene:
                if total > 0:
                    print(last_gene, f"{gene_nucl/total:.1f}", f"{covered/total:.5f}", sep="\t")
                covered = 0
                total = 0
                gene_nucl = 0
                last_gene = cur_gene

            total += span
            gene_nucl += span * read_depth
            if read_depth > 0:
                covered += span


if __name__ == "__main__":
    get_coverage_per_gene(Path(sys.argv[1]))
