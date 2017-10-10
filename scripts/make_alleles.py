#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a bed file into an alleles file.
"""
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse


def bed_to_alleles(mpileup_file, bed_file, allele_file):
    """Convert a bed file into an alleles file.

    Arguments should be writable files or file names
    """
    all_loci = {}
    with open_zipped(bed_file) as fin:
        for line in fin:
            fields = line.strip().split('\t')
            if fields[0] not in all_loci:
                all_loci[fields[0]] = {}
            all_loci[fields[0]][fields[2]] = (fields[5], fields[6])

    with open_zipped(mpileup_file) as fin, open_zipped(allele_file, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            fout.write(
                '{}\t{}\t{}\n'.format(
                    fields[0],
                    fields[1],
                    '\t'.join(
                        all_loci[fields[0]][fields[1]]
                    )
                )
            )
    return 0


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile
    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        else:
            return _bz2.BZ2File(infile, mode)
    return open(infile, mode)


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    # Optional arguments
    parser.add_argument('mpileup',
                        help='The M-pileup file to match')

    # Optional files
    parser.add_argument('bedfile', nargs='?',
                        help="Input file (Default: STDIN)")
    parser.add_argument('allele_file', nargs='?',
                        help="Output file (Default: STDOUT)")

    args = parser.parse_args(argv)

    ifl = args.bedfile if args.bedfile else _sys.stdin
    ofl = args.allele_file if args.allele_file else _sys.stdout

    return bed_to_alleles(args.mpileup, ifl, ofl)

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
