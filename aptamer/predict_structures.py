#!/usr/bin/python

"""
Predicts RNA structures using either ViennaRNA or mfold
(specified in args).
Outputs a fasta with bracket-dot structures.

Overall statistics are printed to standard output.
"""
from __future__ import print_function

import time
import os
import sys
import re
import argparse
from builtins import str
from helpers.functions import *
from helpers.fasta_struct_file import FastaStructFile


def main():
    start_time = time.strftime('%Y-%m-%d %I:%M:%S%p').lower()
    args = parse_arguments()

    in_fname = args.input_file

    fasta_file = FastaFile(in_fname)
    fasta_file.prefix = args.prefix
    fasta_file.suffix = args.suffix
    fasta_file.run_mfold = args.run_mfold
    fasta_file.pass_options = args.pass_options
    fasta_file.vienna_version = args.vienna_version

    # rna_seq_objs = fasta_file

    # output fasta with structure line
    if args.output:
        out_fasta_fname = args.output
    else:
        out_fasta_fname = in_fname + '.struct.fa'

    rna_seq_objs = []
    with open(out_fasta_fname, 'w') as out_fasta_f:
        for node in fasta_file.p_seq_objs():
            out_fasta_f.write(
                '>%s\n%s\n%s\n' % (node.name, node.sequence, node.structure)
            )
            rna_seq_objs.append(node)


    # nodes are now populated. find edges.
    if args.calculate_stats:
        print('Calculating stats...')
        stats = make_aptamer_stats()
        xgmml_obj = XGMML(in_fname)
        find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats)
        print_stats(stats, args)
    else:
        print('Stats calculation disabled')

    print('\n')
    print('Structure fasta written to %s' % (
        out_fasta_fname if (args.calc_structures) else out_xgmml_fname
    ))
    output_stats_tsv(rna_seq_objs, args)
    output_log(args, start_time)
    print()


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
        '-t', '--stats', help=(
            'Tab-separated file to write aptamer structure '
            'statistics to. '
            '(Default: <input_filename>.properties)'
        )
    )
    parser.add_argument(
        '-l', '--log', help=(
            'Log file to which to output command, date, '
            'structure prediction program version, and citation. '
            '(Default: <input_filename>.log)'
        )
    )
    parser.add_argument(
        '-m', '--run_mfold', action='store_true', default=False,
        help=(
            'Run mfold instead of Vienna RNAFold to predict '
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
        '--calculate_stats', action='store_true', default=False,
        help=(
            'If statistics should be calculated. Cannot handle sequences longer than 500'
            '(Default: False)'
        )
    )
    parser.add_argument(
        '--prefix', default='', help=(
            'Sequence to prepend to RNA sequences during '
            'structure prediction. '
            '(Default: --NO PRIMER-- previously used GGGAGGACGAUGCG)'
        )
    )
    parser.add_argument(
        '--suffix', default='', help= (
            'Sequence to append to RNA sequences during '
            'structure prediction. '
            '(Default: --NO PRIMER-- previously used CAGACGACUCGCCCGA)'
        )
    )
    parser.add_argument(
        '-b', '--num_batches', type=int, default=1,
        help=(
            'How many batches to split the input to RNAdistance into'
        )
    )
    parser.add_argument(
        '--pass_options', help=(
            'Quoted string containing options to pass to Vienna or '
            'mfold (depending on whether -m is selected) in lieu '
            'of the following default options. '
            'Vienna default options: "-p -T 30 --noLP --noPS --noGU". '
            'mfold default options: "T=30"'
        )
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        print('\nError: Input file not specified.')
        sys.exit(1)

    # workaround for argparse bug
    # (thinks quoted strings starting with dash are arguments)
    if '--pass_options' in sys.argv:
        pass_options_index = sys.argv.index('--pass_options')
        if pass_options_index < len(sys.argv) - 1:
            next_arg = sys.argv[pass_options_index + 1]
            if (
                (next_arg.startswith('-')) and
                (next_arg not in list(parser._option_string_actions.keys()))
            ):
                sys.argv[pass_options_index + 1] = ' %s' % next_arg

    args = parser.parse_args()

    # validate args
    if args.vienna_version not in ['1', '2', 1, 2]:
        parser.print_help()
        print('Error: ViennaRNA package version not recognized: %s' % (
            args.vienna_version
        ))
        sys.exit(1)

    args.calc_structures = True
    args.edge_type = 'both'
    args.max_edit_dist = 3
    args.max_tree_dist = 3

    return args


def output_stats_tsv(rna_seq_objs, args):
    if args.stats:
        out_fname = args.stats
    else:
        out_fname = args.input_file + '.properties'

    # get header
    categs = list(rna_seq_objs[0].__slots__)
    expand_mfold = False
    if args.run_mfold and ('energy_dict' in categs):
        expand_mfold = True
        energy_dict_keys = list(rna_seq_objs[0].energy_dict.keys())
        categs.remove('energy_dict')
        categs += ['mfold_%s' % z for z in energy_dict_keys]
    # make "name" the first column
    if (categs[0] != 'name') and ('name' in categs):
        name_index = categs.index('name')
        categs[0], categs[name_index] = categs[name_index], categs[0]
    out_list = []
    for x in rna_seq_objs:
        curr_list = []

        for y in categs:
            if expand_mfold and y.replace('mfold_', '') in energy_dict_keys:
                curr_list.append(
                    str(x.energy_dict[y.replace('mfold_', '')])
                )
            elif getattr(x, y) is not None:
                curr_list.append(str(getattr(x, y)))
            else:
                curr_list.append('NA')
        out_list.append(curr_list)

    with open(out_fname, 'w') as out_f:
        out_f.write('%s\n' % '\t'.join(categs))
        for x in out_list:
            out_f.write('%s\n' % '\t'.join(x))
    print('Statistics written to %s' % out_fname)


def output_log(args, start_time):
    if args.log:
        out_fname = args.log
    else:
        out_fname = args.input_file + '.log'
    # get citation
    vienna1_citation = (
        'I.L. Hofacker, W. Fontana, P.F. Stadler, L.S. Bonhoeffer, '
        'M. Tacker, and P. Schuster (1994), "Fast folding and comparison '
        'of RNA secondary structures", Monatshefte fur Chemie 125'
    )
    vienna2_citation = (
        'R. Lorenz, S.H. Bernhart, C. Hoener zu Siederdissen, H. Tafer, C. '
        'Flamm, P.F. Stadler, and I.L. Hofacker (2011), "ViennaRNA Package '
        '2.0", Algorithms for Molecular Biology 6:26'
    )
    mfold_citation = (
        'M. Zuker (1989), "On finding all suboptimal foldings of an RNA '
        'molecule", Science 244:4900'
    )
    if args.run_mfold:
        citation = mfold_citation
    elif str(args.vienna_version) == 1:
        citation = vienna1_citation
    else:
        citation = vienna2_citation

    with open(out_fname, 'w') as out_f:
        out_f.write('Command: %s\n' % ' '.join(sys.argv))
        out_f.write('Start time: %s\n' % start_time)
        if args.run_mfold:
            write_version_str(out_f, 'mfold', 'mfold -v')
        else:
            write_version_str(out_f, 'Vienna RNAFold', 'RNAfold --version')
        out_f.write('Citation: %s\n' % citation)
    print('Log written to %s' % out_fname)


def write_version_str(out_f, program_name, version_command):
    p = subprocess.Popen(
        version_command.split(), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='ascii'
    )
    stdout, stderr = p.communicate()
    if stderr:
        out_f.write(
            'Program: Error getting %s version: %s\n' %
            (program_name, stderr)
        )
    else:
        out_f.write(
            'Program: %s version %s\n' %
            (program_name, re.sub('[^.0-9]', '', stdout))
        )


if __name__ == '__main__':
    sys.exit(main())
