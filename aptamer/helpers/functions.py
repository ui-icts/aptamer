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

import subprocess
import itertools as it
import operator
import textwrap
import numpy
import scipy.stats
from scipy.special import comb as ncr
from Bio import SeqIO
import Levenshtein
import time
import RNA
from multiprocessing import Pool
from helpers.fasta_struct_file import FastaStructFile
from helpers.fasta_file import FastaFile
from helpers.rna_sequence import RNASequence
from helpers.rna_sequence_pair import RNASequencePair
from memory_profiler import profile

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


def find_edges_seed(rna_seq_objs, xgmml_obj, args, stats):
    """Find graph edges using seed algorithm."""
    nodes_copy = list(rna_seq_objs)
    while len(nodes_copy) > 2:
        nodes_copy[0].use_for_comparison = False
        seq_pairs = []
        # go through and find all the matches
        # then remove matches and start again till gone
        for x in range(1, len(nodes_copy)):
            # new nodes_copy[0] each time a node is deleted
            pair = RNASequencePair.build((nodes_copy[0], nodes_copy[x]))
            seq_pairs.append(pair)

        for pair in seq_pairs:
            pair.output(xgmml_obj, args)
            append_pair_stats(stats, pair)

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

def find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats):
    """"Find edges using non-seed algorithm."""

    rna_seq_objs = list(rna_seq_objs)

    num_items = ncr(len(rna_seq_objs), 2)
    batch_size = int(num_items / args.num_batches)
    if args.num_batches > 100:
        batch_size = args.num_batches

    print("{0} items in chunks of {1}".format(num_items, batch_size))

    with Pool() as pool:

        seq_pairs = it.combinations(rna_seq_objs, 2)

        print("STARTING BATCH")

        b_start = time.time()
        result = pool.imap(
                RNASequencePair.build,
                seq_pairs,
                chunksize=batch_size)

        b_end = time.time()
        print("FINISHED BATCH")
        print(b_end - b_start)

        for pair in result:
            pair.output(xgmml_obj, args)
            if args.calculate_stats:
                append_pair_stats(stats, pair)



def find_edges_no_seed_p(rna_seq_objs, xgmml_obj, args, stats):

    #
    # batches_of_pairs look like
    # [(p1,p2,p3,p4),(p5,p6)]
    # and
    # p1 = (s1,s2)
    # so
    # [((s1,s2),(s1,s3),p3,p4),(p5,p6)]
    # and i want to turn it into 
    # [ STRING, STRING,STRING ]
    # But because it would require nested
    # mapping of each of those entries it
    # seems more clear to put that into
    # the rna_distance function itself

    structure = operator.attrgetter('structure')

    rna_seq_objs = list(rna_seq_objs)

    num_items = ncr(len(rna_seq_objs), 2)
    batch_size = int(num_items / args.num_batches)

    if args.num_batches > 1000:
        batch_size = args.num_batches

    print("{0} items in batches of {1}".format(num_items, batch_size))
    just_structures = map(structure, rna_seq_objs)
    pairs_of_structures = it.combinations(just_structures, 2)
    batches_of_pairs = grouper(
            pairs_of_structures,
            batch_size,
            fillvalue=None)

    with Pool() as pool:

        print("STARTING BATCH")
        b_start = time.time()

        batch_results = pool.imap(
                rna_distance,
                batches_of_pairs)


        parsed_results = map(parse_rna_distance_batch, batch_results)
        parsed_results = it.chain.from_iterable(parsed_results)

        seq_pairs = it.combinations(rna_seq_objs, 2)

        pairs_with_results = it.zip_longest(seq_pairs, parsed_results)
        pairs = (RNASequencePair(seq1, seq2, dist) for (seq1, seq2), dist in pairs_with_results)

        b_end = time.time()
        print("FINISHED BATCH")
        print(b_end - b_start)

        for pair in pairs:
            pair.output(xgmml_obj, args)
            if args.calculate_stats:
                append_pair_stats(stats, pair)


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)


def process_seq_pairs(seq_pairs, args):
    """Runs a batch of RNASequence pairs."""

    tree_distance = compute_tree_distances(seq_pairs)

    for i, x in enumerate(seq_pairs):
        x.tree_distance = tree_distance[i]

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

    for stat_pair in it.combinations(stats, 2):
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


def rna_distance(structures):

    # Because the grouper needs to fill in missing
    # values when it uses zip_longest, we need to
    # filter out those values before calling
    # the program

    structures = filter(None, structures)
    structures = it.chain.from_iterable(structures)
    args = '\n'.join(structures)

    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True,
        encoding='ascii'
    )
    stdout_value, stderr_value = sffproc.communicate(args)
    return stdout_value  # f: num1 \n f: num2 \n


def parse_rna_distance_batch(line):
    tree_distance = line.strip('\n')
    return [x.split(':')[1].strip() for x in line.strip('\n').split('\n')]


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




