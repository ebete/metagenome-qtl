#!/usr/bin/env python3

import argparse
import csv
import logging
from pathlib import Path


def get_coverage_per_gene(bed_file):
    logging.info("Reading alignment depths from %s ...", bed_file.name)

    print("orf", "length", "total_mapped", "avg_depth", "coverage", sep="\t")

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
                    print(last_gene, total, gene_nucl, f"{gene_nucl / total:.5f}", f"{covered / total:.5f}", sep="\t")
                covered = 0
                total = 0
                gene_nucl = 0
                last_gene = cur_gene

            total += span
            gene_nucl += span * read_depth
            if read_depth > 0:
                covered += span


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the per-contig coverage and average depths.")
    parser.add_argument("bed_file", metavar="BED", help="Reference alignment depth BED file.")
    args = parser.parse_args()

    logging.basicConfig(format="[%(asctime)s] %(message)s", level=logging.INFO)

    bed_file = Path(args.bed_file)

    try:
        get_coverage_per_gene(bed_file)
    except (BrokenPipeError, KeyboardInterrupt) as ex:
        pass
