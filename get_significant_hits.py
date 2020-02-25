#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

from Bio.Blast import NCBIXML


def get_significant_hits(nxml_file, expect_cutoff):
    logging.info("Getting alignments from %s ...", nxml_file)

    print("query", "has_significant_hit", sep="\t")
    with nxml_file.open("r") as f:
        handle = NCBIXML.parse(f)
        for record in handle:
            query_id = record.query_id
            has_significant_hit = False
            for alignment in record.alignments:
                # id_acc = alignment.hit_id
                e_value = alignment.hsps[0].expect
                if e_value > expect_cutoff:
                    continue
                has_significant_hit = True
            print(query_id, "found" if has_significant_hit else "not_found", sep="\t")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")
    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_nxml", help="NCBI-XML BLAST results file.", metavar="NXML")
    parser.add_argument("-e", "--expect", help="Minimum e-value.", metavar="E", type=float, default=1e-3)
    args = parser.parse_args()

    xml_path = Path(args.input_nxml).expanduser().resolve().absolute()
    expect_cutoff = args.expect

    get_significant_hits(xml_path, expect_cutoff)
