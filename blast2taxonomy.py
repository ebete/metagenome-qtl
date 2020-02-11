#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

from Bio.Blast import NCBIXML


def read_nxml(xml_path):
    print("query", "hit_id", "identity", "align_length", "expect", sep="\t")
    with xml_path.open("r") as f:
        handle = NCBIXML.parse(f)
        for record in handle:
            if len(record.alignments) == 0:
                print(record.query, "NA", "NA", "NA", "NA", sep="\t")

            for alignment in record.alignments:
                id_acc = alignment.title
                e_value = alignment.hsps[0].expect
                aln_length = alignment.hsps[0].align_length
                identity = alignment.hsps[0].identities
                print(record.query, id_acc, identity, aln_length, e_value, sep="\t")

                print("#query: ", alignment.hsps[0].query)
                print("#       ", alignment.hsps[0].match)
                dot_format = ""
                for i, m in enumerate(alignment.hsps[0].match):
                    dot_format += "." if m == "|" else alignment.hsps[0].sbjct[i]
                print("#target:", dot_format)
                break
            print("#")


if __name__ == "__main__":
    logging.basicConfig(format="[%(asctime)s] %(message)s", level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("nxml_file", metavar="NXML", help="NCBI-XML BLAST results file.")
    args = parser.parse_args()

    read_nxml(Path(args.nxml_file))
