#!/usr/bin/env python3

"""
Copyright (c) 2020 Thom Griffioen <t.griffioen@nioo.knaw.nl>
MIT License

Created: 2020-01-31
"""

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path

import pysam
from Bio import SeqIO


def get_mapped_reads(samfile, mapq_thresh=0):
    """
    Retrieves a set of all aligned read IDs.

    :type samfile: pathlib.Path
    :param samfile: Input SAM file.

    :type mapq_thresh: int
    :param mapq_thresh: Minimum MAPQ to be considered mapped.

    :rtype: dict[str, set[str]]
    :return: Set containing all aligned read IDs.
    """
    read_ids = {
        "fwd": set(),
        "rev": set()
    }

    logging.info(f"Reading SAM file {samfile} ...")
    with pysam.AlignmentFile(samfile, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                logging.debug(f"Skipping {read.qname} (unaligned)")
                continue
            if read.mapq < mapq_thresh:
                logging.debug(f"Skipping {read.qname} (low MQ)")
                continue

            read_ids["fwd" if read.is_read1 else "rev"].add(read.qname)

    logging.info(f"Found {len(read_ids['fwd'])} forward and {len(read_ids['rev'])} reverse reads mapped, "
                 f"{len(read_ids['fwd'] & read_ids['rev'])} paired.")

    return read_ids


def extract_mapped_read_records(fastq_fwd, fastq_rev, read_ids, paired_out, unpaired_out, unaligned_out):
    """
    Split FASTQ files into a FASTQ with the aligned and one with the unaligned reads.

    :param fastq_fwd: Forward FASTQ file used for creating the SAM file.
    :param fastq_rev: Reverse FASTQ file used for creating the SAM file.
    :param read_ids: Read IDs that passed the MAPQ threshold.
    :param paired_out: FASTQ file to write paired aligned reads to.
    :param unpaired_out: FASTQ file to write singleton aligned reads to.
    :param unaligned_out: FASTQ file to write unaligned reads to.
    """
    logging.info("Getting FASTQ records that have mapped fragments ...")

    fq_format = "fastq"
    mappable = 0
    records = 0
    with gzip.open(fastq_fwd, "rt") as fwd_handle, gzip.open(fastq_rev, "rt") as rev_handle:
        for (fwd_record, rev_record) in zip(SeqIO.parse(fwd_handle, fq_format), SeqIO.parse(rev_handle, fq_format)):
            records += 2
            if fwd_record.id in read_ids["fwd"] and rev_record.id in read_ids["rev"]:
                mappable += 2
                # write interleaved
                SeqIO.write(fwd_record, paired_out, fq_format)
                SeqIO.write(rev_record, paired_out, fq_format)
                continue
            if fwd_record.id in read_ids["fwd"]:
                mappable += 1
                SeqIO.write(fwd_record, unpaired_out, fq_format)
            else:
                SeqIO.write(fwd_record, unaligned_out, fq_format)
            if rev_record.id in read_ids["rev"]:
                mappable += 1
                SeqIO.write(rev_record, unpaired_out, fq_format)
            else:
                SeqIO.write(rev_record, unaligned_out, fq_format)
    logging.info(f"{mappable:d} records extracted from {records:d} total ({mappable / records * 100:.1f}%).")


def init_logger(log_level=logging.INFO):
    """
    Initialise the logger.

    :type log_level: int
    :param log_level: Set the minimum logging level.
    """
    logging.basicConfig(level=log_level, format="[%(asctime)s] %(message)s")


if __name__ == '__main__':
    init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="SAM/BAM file to extract mapped fragments from.", metavar="SAM")
    parser.add_argument("input_fastq_fwd", help="Forward FASTQ file with all the query sequences used in generating "
                                                "the SAM file (gzipped).", metavar="FWD_FQ")
    parser.add_argument("input_fastq_rev", help="Reverse FASTQ file with all the query sequences used in generating "
                                                "the SAM file (gzipped).", metavar="REV_FQ")
    parser.add_argument("-u", "--unmapped", help="Output all unmapped sequences to this file instead of stdout.",
                        metavar="UNMAPPED", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("-p", "--paired", help="Optionally output all paired mappable reads to this file.",
                        metavar="PAIRED", type=argparse.FileType("w"), default=os.devnull)
    parser.add_argument("-s", "--singletons", help="Optionally output all unmapped sequences to this file.",
                        metavar="SINGLE", type=argparse.FileType("w"), default=os.devnull)
    parser.add_argument("-q", "--mapq", help="Set the minimum MAPQ score to be considered mapped.",
                        metavar="MQ", type=int, default=0)
    args = parser.parse_args()

    sam_file = Path(args.input_sam).expanduser().resolve().absolute()
    fq_fwd = Path(args.input_fastq_fwd).expanduser().resolve().absolute()
    fq_rev = Path(args.input_fastq_rev).expanduser().resolve().absolute()
    unmapped_fastq = args.unmapped
    paired_fastq = args.paired
    unpaired_fastq = args.singletons
    min_mapq = args.mapq

    try:
        mapped_read_ids = get_mapped_reads(sam_file)
        extract_mapped_read_records(fq_fwd, fq_rev, mapped_read_ids, paired_fastq, unpaired_fastq, unmapped_fastq)
    finally:
        args.paired.close()
        args.singletons.close()
        args.unmapped.close()

    logging.shutdown()
