"""
Functions for predicting RNA structures, and creating graph from a set
of RNA sequences with known structures.

If RNA structures are not already known, outputs a fasta structure file
using structures determined via either ViennaRNA or mfold.

Otherwise, outputs an xgmml graph file in which vertices are RNA 
sequences. Edges are created between vertices that have a small enough
tree distance or edit distance between them.

Overall statistics are printed to standard output.
"""
from __future__ import print_function

from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import os
import sys
import re
import subprocess
import itertools
import textwrap
import argparse
import numpy
import scipy.stats
from Bio import SeqIO
import Levenshtein
import time
import pdb
import RNA
import concurrent.futures
from multiprocessing import Pool
from helpers.fasta_struct_file import FastaStructFile
from helpers.fasta_file import FastaFile
from helpers.rna_sequence import RNASequence
from helpers.rna_sequence_pair import RNASequencePair


class XGMML(object):
    """Graph. XGMML file."""
    def __init__(self, name):
        self.name = name
        self.out_str = ''
        self.nodes = {}

        # [source ID, target ID, (attribute1 name, attribute1 value),
        # (attribute2 name, attribute2 value), ...]
        self.edges = []

    def output_att(self, type_, name, label, value):
        self.out_str += (
            '<att type="%s" name="%s" label="%s" value="%s"/>\n' % (
                type_, name, label, value
            )
        )

    def output(self, args):
        self.out_str = textwrap.dedent(
            """
            <?xml version="1.0"?>
            <graph directed="1" id="5" label="%s"
            xmlns:xsi="http:///www.w3.org/2001/XMLSchema-instance"
            xmlns:ns1="http://www.w3.org/1999/xlink"
            xmlns:dc="http://purl.org/dc/elements/1.1/"
            xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
            xmlns="http://www.cs.rpi.edu/XGMML">
            """ % self.name
        ).lstrip()
        for n in self.nodes:
            self.out_str += '<node id="%s" label="%s" weight="%s">\n' % (
                n, self.nodes[n].name, self.nodes[n].cluster_size
            )
            if not self.nodes[n] == None:
                self.output_att(
                    'integer', 'size', 'Size', self.nodes[n].cluster_size
                )
            self.output_att(
                'string', 'structure', 'Structure', self.nodes[n].structure
            )
            self.output_att(
                'string', 'sequence', 'Sequence', self.nodes[n].sequence
            )
            self.output_att(
                'real', 'energy', 'Energy', self.nodes[n].free_energy
            )
            self.output_att(
                'real', 'ensembleFreeEnergy', 'ensemble Free Energy',
                self.nodes[n].ensemble_free_energy
            )
            self.output_att(
                'real', 'ensembleProbability', 'ensemble Probability',
                self.nodes[n].ensemble_probability
            )
            self.output_att(
                'real', 'ensembleDiversity', 'ensemble Diversity',
                self.nodes[n].ensemble_diversity
            )
            self.out_str += '</node>\n'

        for e in self.edges:
            source, target = e[:2]
            atts = e[2:]
            self.out_str += (
                '<edge source="%s" target="%s" label="%s to %s">\n' % (
                    source, target, source, target
                )
            )
            for att_type, att_name, att_label, att_value in atts:
                self.output_att(att_type, att_name, att_label, att_value)
            self.out_str += '</edge>\n'
        self.out_str += '</graph>\n'
        return self.out_str






def rna_distance(structures):

    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True,
        encoding='ascii'
    )
    stdout_value, stderr_value = sffproc.communicate(structures)
    return stdout_value  # f: num1 \n f: num2 \n


def compute_tree_distances(seq_pairs):
    # Collect all the structures into a single string
    # so that it can be piped to RNADistance
    seq_to_tree = []
    for x in seq_pairs:
        seq_to_tree.append(x.sequence1.structure)
        seq_to_tree.append(x.sequence2.structure)
    string_of_sequences = '\n'.join(seq_to_tree)

    tree_distance = rna_distance(string_of_sequences)
    tree_distance = tree_distance.strip('\n').split('\n')  # take off last lr

    assert len(tree_distance) == len(seq_pairs), (
        'Error length of tree distance %s does not match length of seq_pairs '
        '%s and should -- check installation of RNAdistance'
        % (len(tree_distance), len(seq_pairs))
    )

    return [td.split(' ')[1] for td in tree_distance]


def process_seq_pairs(seq_pairs, args):
    """Runs a batch of RNASequence pairs."""

    tree_distance = compute_tree_distances(seq_pairs)

    for i, x in enumerate(seq_pairs):
        x.tree_distance = tree_distance[i]


def get_mfold_stats(det_filename):
    # fixed values within line
    valid_pattern = {
        0: 'dG', 1: '=', 3: 'dH', 4: '=', 6: 'dS', 7: '=', 9: 'Tm', 10: '='
    }
    energy_stats = {}
    for key, value in valid_pattern.items():
        if value != '=':
            energy_stats[value] = None
    with open(det_filename, 'r') as f:
        for i, row in enumerate(f):
            if i == 5: #6th line
                row = row.split()
                test_pattern = [
                    row[z] ==  valid_pattern[z] for z in valid_pattern
                ]
                assert False not in test_pattern, (
                    'mfold file *.txt.det does not match the expected format'
                )
                for x in valid_pattern:
                    if row[x] != '=':
                        energy_stats[row[x]] = row[x + 2] 
    return energy_stats


def convert_ct_to_bracket_dot(ct_filename):
    bracket_dot = ''
    with open(ct_filename, 'r') as f:
        for row in f:
            row = row.split()
            # used to grab energy but not needed except to skip first line
            if '=' in row and len(row) < 6:  # first row
                pass
            elif row[4] == '0':
                bracket_dot += '.'
            elif int(row[0]) < int(row[4]):
                bracket_dot += '('
            elif int(row[0]) > int(row[4]):
                bracket_dot += ')'
    return '%s' % (bracket_dot) if (bracket_dot != '') else None


def find_edges_seed(rna_seq_objs, xgmml_obj, args, stats):
    """Find graph edges using seed algorithm."""
    nodes_copy = rna_seq_objs
    while len(nodes_copy) > 2:
        nodes_copy[0].use_for_comparison = False
        seq_pairs = []
         # go through and find all the matches
         # then remove matches and start again till gone
        for x in range(1, len(nodes_copy)):
            # new nodes_copy[0] each time a node is deleted
            pair = RNASequencePair( nodes_copy[0], nodes_copy[x])
            pair.energy_delta = abs(
                pair.sequence1.free_energy - pair.sequence2.free_energy
            )
            pair.edit_distance = Levenshtein.distance(
                pair.sequence1.sequence, pair.sequence2.sequence
            )
            seq_pairs.append(pair)
        process_seq_pairs(seq_pairs, args)
        for x in seq_pairs:
            if x.is_valid_edge:
                x.sequence2.use_for_comparison = False
        # delete nodes which already belong to an edge
        new_nodes_copy = [
            z for z in nodes_copy if z.use_for_comparison
        ]
        print('Number of RNA sequences reduced from %d to %d ' % (
            len(nodes_copy), len(new_nodes_copy)
        ))
        nodes_copy = new_nodes_copy

def pairwise_combine(rna_seq_objs):
    """
    Combines the RNASequence objects into pairs.
    A B C D -> [A B], [A C], [A D], [B, C] ...
    """
    for i in range(0, len(rna_seq_objs)):
        for j in range(i + 1, len(rna_seq_objs)):
            pair = RNASequencePair(
                    rna_seq_objs[i], rna_seq_objs[j]
                    )
            yield pair


def grouper(n, iterable):
    """
    >>> list(grouper(3, 'ABCDEFG'))
    [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    """
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])


def find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats):
    """"Find edges using non-seed algorithm."""

    def my_callback(pairs):
        for pair in pairs:
            pair.output(xgmml_obj, args)
            append_pair_stats(stats, pair)

    pool = Pool()

    seq_pairs = itertools.combinations(rna_seq_objs, 2)
    result = pool.starmap_async(
            RNASequencePair,
            seq_pairs,
            callback=my_callback)
    pool.close()
    pool.join()


def make_aptamer_stats():
    return {'energy_delta': [], 'edit_distance': [], 'tree_distance': []}


def append_pair_stats(stats, pair):
    """
    Appends the stat elements from the pair
    to the lists in the stats
    """
    stats['energy_delta'].append(pair.energy_delta)
    stats['edit_distance'].append(pair.edit_distance)
    try:
        stats['tree_distance'].append(float(pair.tree_distance))
    except ValueError:
        stats['tree_distance'].append(None)


def print_stats(stats, args):
    print('\nOverall Statistics:')
    print('--------------------')
    for stat in stats:
        stat_label = stat.capitalize().replace('_', ' ')
        if any(stats[stat]):
            print('%s mean: %.3g' % (stat_label, numpy.mean(stats[stat])))
            print('%s SD: %.3g' % (stat_label, numpy.std(stats[stat])))
            print('%s SEM: %.3g' % (stat_label, scipy.stats.sem(stats[stat])))
        else:
            print('%s mean: N/A' % (stat_label))
            print('%s SD: N/A' % (stat_label))
            print('%s SEM: N/A' % (stat_label))
        print()

    for stat_pair in itertools.combinations(stats, 2):
        print('Pearson correlation(%s, %s):' % (
            stat_pair[0], stat_pair[1]
        ), end=' ')
        if args.calc_structures:
            pearson_coeff, p_value = scipy.stats.pearsonr(
                stats[stat_pair[0]], stats[stat_pair[1]]
            )
            print('%.3g, p=%.3g' % (pearson_coeff, p_value))
        else:
            print('N/A')