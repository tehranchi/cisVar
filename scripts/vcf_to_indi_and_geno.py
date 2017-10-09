#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parse VCF files into genotypes and individuals files.

EX:
    chr   position    ref alt NA18486 NA18489 ...
    chr1  10583   g   a   0   1 ...
    chr1  10611   C   G   0   0 ...
    chr1  13302   c   t   1   0 …

and:
    NA18486
    NA18489
"""
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse

from tqdm import tqdm


EXPECTED_HEADER = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
]


def parse_file(vcfs, geno, indi, logfile=_sys.stderr):
    """File parser."""
    if not isinstance(vcfs, list):
        raise ValueError('vcfs must be list')
    # Tracking variables
    header = None
    idone = False  # Individual order parsing
    _sys.stdout.write('Looping through {} VCF files\n'.format(len(vcfs)))
    with _open_zipped(geno, 'w') as gn:
        for vcf in tqdm(vcfs, unit='file'):
            with _open_zipped(vcf) as fin:
                hdone = False  # Header parsing
                fdone = False  # Format string parsing
                for line in fin:
                    if line.startswith('#'):
                        header = line
                        continue
                    if not hdone:
                        header = header.strip().split('\t')
                        assert header[:9] == EXPECTED_HEADER
                        inds = header[9:]
                        if not idone:
                            orig_inds = inds
                            # Write individual file
                            with _open_zipped(indi, 'w') as ind:
                                ind.write('\n'.join(inds) + '\n')
                            # Write header
                            gn.write('chr\tposition\tref\talt\t{}\n'.format(
                                '\t'.join(inds)
                            ))
                            idone = True
                        if not inds == orig_inds:
                            raise ValueError(
                                ("Individual order is different in {fl}!\n"
                                "Orig: {orig}\nCurr: {cur}").format(
                                    fl=vcf, orig=orig_inds, cur=inds
                                )
                            )
                        hdone = True
                    f = line.strip().split('\t')
                    if not fdone:
                        # Get the genotype info location
                        try:
                            gti = f[8].split(':').index('GT')
                        except ValueError:
                            raise ValueError(
                                'GT column missing from VCF format string. In need '
                                'the genotype info location to know where to look '
                                'for the genotype (looks like 0|1)'
                            )
                    # Write core records
                    outstr = (
                        '{chr}\t{pos}\t{ref}\t{alt}'.format(
                            chr=f[0], pos=f[1], ref=f[3], alt=f[4]
                        )
                    )
                    # Write genotypes
                    badgt = False
                    for ind in f[9:]:
                        gt = ind.split(':')[gti]
                        if not '|' in gt:
                            raise ValueError(
                                'Expected genotype to be like 0|0, is {}'.format(gt) +
                                '\n' + line
                            )
                        if gt == '0|0':
                            gt = '0'
                        elif gt in ['0|1', '1|0']:
                            gt = '1'
                        elif gt == '1|1':
                            gt = '2'
                        else:
                            logfile.write(
                                'Invalid genotype {}: {}'.format(gt, line)
                            )
                            gt = 'NA'
                            badgt = True
                            break
                        outstr += ('\t{}'.format(gt))
                    if badgt:
                        continue
                    gn.write(outstr + '\n')
    return


def _open_zipped(infile, mode='r'):
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

    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('genotype_file', help="Genotype file")
    parser.add_argument('individual_file', help="Genotype file")
    parser.add_argument('vcf_files', nargs='+', help="VCF Files")

    parser.add_argument('-l', '--logfile', default=_sys.stderr,
                        type=_argparse.FileType('a'),
                        help="Log file (Default: STDERR)")

    args = parser.parse_args(argv)

    parse_file(
        args.vcf_files, args.genotype_file,
        args.individual_file, args.logfile
    )


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
