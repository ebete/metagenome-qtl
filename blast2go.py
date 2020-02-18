#!/usr/bin/env python3

import argparse
import logging
import lzma
import pickle
import sys
from pathlib import Path

from Bio.Blast import NCBIXML


def get_id_lookup_table(serialised_db):
    logging.info("Loading ID lookup file %s", serialised_db)
    with lzma.open(serialised_db, "r") as f:
        return pickle.load(f)


def read_nxml(xml_path, lookup_db_path):
    lookup_db = get_id_lookup_table(lookup_db_path)

    with xml_path.open("r") as f:
        handle = NCBIXML.parse(f)
        for record in handle:
            for alignment in record.alignments:
                id_acc = alignment.hit_id
                e_value = alignment.hsps[0].expect
                print(record.query, id_acc, "; ".join(lookup_db.get(id_acc, ["None"])), e_value, sep="\t")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_lookup", help="ID conversion file generated with "
            "make_uniprot_idmapping_db.py", metavar="LOOKUP_DB")
    parser.add_argument("nxml_file", help="NCBI-XML BLAST results file.", metavar="NXML")
    args = parser.parse_args()

    input_lookup = Path(args.input_lookup)
    nxml_file = Path(args.nxml_file)

    read_nxml(nxml_file, input_lookup)
