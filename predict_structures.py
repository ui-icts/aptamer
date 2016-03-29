#!/usr/bin/python

"""
Predicts RNA structures using either ViennaRNA or mfold
(specified in args).
Outputs a fasta with bracket-dot structures.

Overall statistics are printed to standard output.
"""

import os
import sys
import re
import subprocess
import itertools
import textwrap
import argparse
import numpy
import scipy.stats
import Levenshtein

from aptamer_functions import *


def main():
    args = parse_arguments()
    cluster_size_re = re.compile('SIZE=(\d+)')
    in_fname = args.input_file
    in_fh = open(in_fname)
    stats = {'energy_delta':[], 'edit_distance':[], 'tree_distance':[]}
    rna_seq_objs = []  # list of RNASequence objects

    process_fasta(in_fh, args, cluster_size_re, rna_seq_objs)
    # output fasta with structure line
    if args.output:
        out_fasta_fname = args.output
    else:
        out_fasta_fname = in_fname + '.struct.fa'
    with open(out_fasta_fname, 'w') as out_fasta_f:
        for node in rna_seq_objs:
            out_fasta_f.write(
                '>%s\n%s\n%s\n' % (node.name, node.sequence, node.structure)
            )
    process_struct_fasta(in_fh, args, cluster_size_re, rna_seq_objs)
    xgmml_obj = XGMML(in_fname)

    # nodes are now populated. find edges.
    print 'Calculating stats...'
    find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats)
    print_stats(stats, args)

    print '\n'
    print 'Output written to %s.' % (
        out_fasta_fname if (args.calc_structures) else out_xgmml_fname
    )
    output_stats_tsv(rna_seq_objs, args)
    print
    in_fh.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            '%s: Predict RNA structures using either ViennaRNA or mfold '
            '(specified in args). '
            'Output a fasta with bracket-dot structures.' % (
                os.path.basename(sys.argv[0])
            )
        )
    )
    parser.add_argument(
        'input_file', help=(
            'Input fasta file containing RNA sequences. '
            'Must be valid fasta without a bracket-dot structure line.'
        )
    )
    parser.add_argument(
        '-o', '--output', help=(
            'Specify path of output fasta file with structure line. '
            '(Default: <input_filename>.struct.fa)'
        )
    )
    parser.add_argument(
        '-m', '--run_mfold', action='store_true', default=False,
        help=(
            'Run mfold instead of to Vienna RNAFold to predict '
            'RNA structures. Note: some graph attributes will '
            'be omitted.'
        )
    )
    parser.add_argument(
        '-v', '--vienna_version', type=int, default=2,
        help=(
            'ViennaRNA package version. Specify "1" or "2". '
            '(Default: 2)'
        )
    )
    parser.add_argument(
        '-p', '--prefix', default='', help=(
            'Sequence to prepend to RNA sequences during '
            'structure prediction. '
            '(Default: --NO PRIMER-- previously used GGGAGGACGAUGCG)'
        )
    )
    parser.add_argument(
        '-s', '--suffix', default='', help= (
            'Sequence to append to RNA sequences during '
            'structure prediction. '
            '(Default: --NO PRIMER-- previously used CAGACGACUCGCCCGA)'
        )
    )
    parser.add_argument(
        '-t', '--stats_file', help=(
            'Tab-separated file to write aptamer structure '
            'statistics to. '
            '(Default: <input_filename>.stats.tsv)'
        )
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        print '\nError: Input file not specified.'
        sys.exit(1)

    args = parser.parse_args()

    # validate args
    if args.vienna_version not in ['1', '2', 1, 2]:
        parser.print_help()
        print 'Error: ViennaRNA package version not recognized: %s' % (
            args.vienna_version
        )
        sys.exit(1)

    args.calc_structures = True
    args.edge_type = 'both'
    args.max_edit_dist = 3
    args.max_tree_dist = 3

    return args


def output_stats_tsv(rna_seq_objs, args):
    if args.stats_file:
        out_fname = args.stats_file
    else:
        out_fname = args.input_file + '.stats.tsv'
    categs = rna_seq_objs[0].__dict__.keys()

    # make "name" the first column
    if (categs[0] != 'name') and ('name' in categs):
        name_index = categs.index('name')
        categs[0], categs[name_index] = categs[name_index], categs[0]

    out_list = []
    for x in rna_seq_objs:
        curr_list = []
        d = x.__dict__
        for y in categs:
            if y in d and d[y]:
                curr_list.append(str(d[y]))
            else:
                curr_list.append('NA')
        out_list.append(curr_list)

    #if os.path.exists(out_fname):
    #    print (
    #        'Warning: Overwriting existing stats file "%s".' % out_fname
    #    )
    with open(out_fname, 'w') as out_f:
        out_f.write('%s\n' % '\t'.join(categs))
        for x in out_list:
            out_f.write('%s\n' % '\t'.join(x))
    print 'Statistics written to %s.' % out_fname


if __name__ == '__main__':
    sys.exit(main())
