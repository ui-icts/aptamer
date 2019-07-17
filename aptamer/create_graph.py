#!/usr/bin/python

"""
Creates a graph from a set of RNA sequences with already-predicted
structures.

Outputs an xgmml graph file in which vertices are RNA 
sequences. Edges are created between vertices that have a small enough
tree distance or edit distance between them. (Specified in args.)
"""
from __future__ import print_function

from helpers.functions import *
from helpers.fasta_struct_file import FastaStructFile
import argparse
import time
import os
import sys

def main():
    start_time = time.strftime('%Y-%m-%d %I:%M:%S%p').lower()
    args = parse_arguments()
    in_fname = args.input_file

    if args.save_combinations:
        with open(in_fname, 'r') as file_handle:
            FastaStructFile(file_handle).write_combinations("combinations.txt")

    in_fh = open(in_fname, 'r')
    struct_file = FastaStructFile(in_fh)
    rna_seq_objs = struct_file.rna_seq_objs()

    xgmml_obj = XGMML(in_fname)

    stats = make_aptamer_stats()
    # nodes are now populated. find edges.
    if args.seed:
        print('Finding edges (seed)...')
        find_edges_seed(rna_seq_objs, xgmml_obj, args, stats)
    else:
        print('Finding edges (no seed)...')
        if args.spawn:
            find_edges_no_seed_p(rna_seq_objs, xgmml_obj, args, stats)
        else:
            find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats)


    # output xgmml file
    if args.output:
        out_xgmml_fname = args.output
    else:
        out_xgmml_fname = in_fname + '.xgmml'
    with open(out_xgmml_fname, 'w') as out_xgmml_f:
        out_xgmml_f.write(xgmml_obj.output(args))

    if args.calculate_stats:
        print_stats(stats, args)

    if args.calculate_stats:
        print('\n\nOutput written to %s' % (
            out_fasta_fname if (args.calc_structures) else out_xgmml_fname
        ))

    output_log(args, start_time)
    print()
    in_fh.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            '%s: Create a graph from a set of RNA sequences with '
            'already-predicted structures.' % os.path.basename(sys.argv[0])
        )
    )
    parser.add_argument(
        'input_file', help=(
            'Input fasta containing RNA sequences. '
            'Must contain a bracket-dot structure line.'
        )
    )
    parser.add_argument(
        '-o', '--output', help=(
            'Specify path of output xgmml graph file. '
            '(Default: <input_filename>.xgmml)'
        )
    )
    parser.add_argument(
        '-l', '--log', help=(
            'Log file to which to output command and date. '
            '(Default: <input_filename>.log)'
        )
    )
    parser.add_argument(
        '-t', '--edge_type', default='both', help=(
            'Whether to create edges in output graph according to '
            'edit distance or tree distance.'
            '(Specify "edit", "tree", or "both".) (Default: both)'
        )
    )
    parser.add_argument(
        '-e', '--max_edit_dist', type=int, default=3,
        help=(
            'Maximum edit distance allowed for an edge to be created '
            'in output graph. '
            'Assumes --edge_type is "edit" or "both". (Default: 3)'
        )
    )
    parser.add_argument(
        '-d', '--max_tree_dist', type=int, default=3,
        help=(
            'Maximum tree distance allowed for an edge to be created '
            'in output graph. '
            'Assumes --edge_type is "tree" or "both". (Default: 3)'
        )
    )
    parser.add_argument(
        '-b', '--num_batches', type=int, default=1,
        help=(
            'How many batches to split the input to RNAdistance into'
        )
    )
    parser.add_argument(
        '--seed', action='store_true', default=False,
        help='Use seed sequence algorithm to find graph edges.'
    )
    parser.add_argument(
        '--save_combinations', action='store_true', default=False,
        help='Save the combinations used to a file for inspection.'
    )
    parser.add_argument(
        '--spawn', action='store_true', default=False,
        help='Run RNAdistance by spawning and external exe'
    )
    parser.add_argument(
        '--calculate_stats', action='store_true', default=False,
        help=(
            'If statistics should be calculated. Cannot handle sequences longer than 500'
            '(Default: False)'
        )
    )

    if len(sys.argv) <= 1:
        parser.print_help()
        print('\nError: Input file not specified.')
        sys.exit(1)

    args = parser.parse_args()

    # validate args
    args.edge_type = args.edge_type.lower()
    if args.edge_type not in ['edit', 'tree', 'both']:
        parser.print_help()
        print('Error: Edge type option not recognized: %s' % (
            args.edge_type
        ))
        sys.exit(1)

    args.calc_structures = False

    return args


def output_log(args, start_time):
    if args.log:
        out_fname = args.log
    else:
        out_fname = args.input_file + '.log'

    with open(out_fname, 'w') as out_f:
        out_f.write('Command: %s\n' % ' '.join(sys.argv))
        out_f.write('Start time: %s\n' % start_time)

    if args.calculate_stats:
        print('Log written to %s' % out_fname)


if __name__ == '__main__':
    sys.exit(main())
