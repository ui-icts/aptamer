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

class RNASequence:
    """Graph node."""
    def __init__(self, name, cluster_size, seq):
        self.name = name  # integer ID
        self.cluster_size = int(cluster_size)
        self.sequence = seq
        self.structure = None
        self.free_energy = None
        self.energy_dict = {}
        self.ensemble_free_energy = None
        self.ensemble_probability = None
        self.ensemble_diversity = None
        self.use_for_comparison = True  # used only with seed option

    def full_output(self):
        attrs = vars(self)
        print ','.join('%s:%s' % item for item in attrs.items())

    def output(self):
        print '>%s  ' % (self.name)
        print self.sequence
        print self.structure

    def __str__(self):
        return '\n'.join(
            ['%s : %s' % (z, self.__dict__[z]) for z in self.__dict__]
        )


class RNASequencePair:
    """Graph edge. Pair of RNASequence objects."""
    def __init__(self, seq1, seq2):
        self.sequence1 = seq1
        self.sequence2 = seq2
        self.energy_delta = None
        self.edit_distance = None
        self.tree_distance = None
        self.is_valid_edge = False  # used only with seed option

    def __str__(self):
        return '%s\n---\n%s' % (str(self.sequence1), str(self.sequence2))

    def output(self, xgmml, args):
        # if the xgmml data structure does not have this node, add it
        if self.sequence1.name not in xgmml.nodes:  
            xgmml.nodes[self.sequence1.name] = self.sequence1
        if self.sequence2.name not in xgmml.nodes:
            xgmml.nodes[self.sequence2.name] = self.sequence2

        # make edge between nodes that are similar enough
        # in terms of edit or tree distance
        interaction = []
        if args.edge_type in ['edit', 'both']:
            if (
                self.edit_distance and
                (int(self.edit_distance) <= args.max_edit_dist)
            ):
                interaction.append('edit distance')
                self.is_valid_edge = True
        if args.edge_type in ['tree', 'both']:
            if (
                self.tree_distance and
                (int(self.tree_distance) <= args.max_tree_dist)
            ):
                interaction.append('tree distance')
                self.is_valid_edge = True

        if not self.is_valid_edge:
            return
        if len(interaction) > 1:
            interaction = 'both'
        else:
            interaction = interaction[0]

        xgmml.edges.append([
            self.sequence1.name, self.sequence2.name,
            ('string', 'interaction', 'interaction', interaction),
            ('integer', 'editDistance', 'edit distance', self.edit_distance),
            ('integer', 'treeDistance', 'tree distance', self.tree_distance)
        ])


class XGMML:
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


def run_rnafold(seq, args):
    print '##################'
    print 'Running RNAFold...'
    print '##################'
    cmd = None
    if args.pass_options is not None:
        cmd = ['RNAfold %s' % args.pass_options]
    elif args.vienna_version == 1:
        cmd = ['RNAfold -p -T 30 -noLP -noPS -noGU']
    elif args.vienna_version == 2:
        cmd = ['RNAfold -p -T 30 --noLP --noPS --noGU']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True
    )
    print 'Command:', cmd[0]
    stdout_value, stderr_value = sffproc.communicate(seq)
    return stdout_value


def run_mfold(seq, args):
    print '##################'
    print 'Running mfold...'
    print '##################'
    if not os.path.exists('mfold_out'):
        os.mkdir('mfold_out')
    os.chdir('mfold_out')
    temp_filename = 'mfold_temp.txt'
    with open(temp_filename, 'w') as f:
        f.write('%s\n' % seq)
    if args.pass_options is not None:
        cmd_string = 'mfold SEQ=%s %s' % (temp_filename, args.pass_options)
    else:
        cmd_string = 'mfold SEQ=%s T=30' % temp_filename
    ret = subprocess.call(
        cmd_string, stderr=subprocess.STDOUT, shell=True
    )
    if ret != 0:
        print (
            'Error when running mfold. Return code %d. '
            'See mfold log file for details.' % ret
        )
        sys.exit(ret)
    print
    structure = convert_ct_to_bracket_dot('%s.ct' % temp_filename)
    energy_stats = get_mfold_stats('%s.det' % temp_filename)
    os.chdir('..')
    return structure, energy_stats


def rna_distance(structures):
    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True
    )
    stdout_value, stderr_value = sffproc.communicate(structures)
    return stdout_value #f: num1 \n f: num2 \n


def process_seq_pairs(seq_pairs, xgmml, args, stats):
    """Runs a batch of RNASequence pairs."""
    #print '\n\n'.join([str(z) for z in seq_pairs])
    seq_to_tree = []
    append_s = time.clock()
    for x in seq_pairs:
        seq_to_tree.append(x.sequence1.structure)
        seq_to_tree.append(x.sequence2.structure)
    append_e = time.clock()

    tree_distance = rna_distance('\n'.join(seq_to_tree))
    tree_distance = tree_distance.strip('\n').split('\n')  # take off last lr
    distance_e = time.clock()

    assert len(tree_distance) == len(seq_pairs), (
        'Error length of tree distance %s does not match length of seq_pairs '
        '%s and should -- check installation of RNAdistance'
        % (len(tree_distance), len(seq_pairs))
    )

    for i, x in enumerate(seq_pairs):
        # TODO: I think seq_pairs[i] == x here
        # should just use one or the other
        seq_pairs[i].tree_distance = tree_distance[i].split(' ')[1]
        x.output(xgmml, args)
        stats['energy_delta'].append(x.energy_delta)
        stats['edit_distance'].append(x.edit_distance)
        try:
            stats['tree_distance'].append(float(x.tree_distance))
        except ValueError:
            stats['tree_distance'].append(None)

    output_e = time.clock()
    append_time = append_e - append_s
    distance_time = distance_e - append_e
    output_time = output_e - distance_e
    print 'APPEND:   {}'.format(append_time)
    print 'DISTANCE: {}'.format(distance_time)
    print 'OUTPUT:   {}'.format(output_time)

def get_mfold_stats(det_filename):
    #fixed values within line
    valid_pattern = {
        0:'dG', 1:'=', 3:'dH', 4:'=', 6:'dS', 7:'=', 9:'Tm', 10:'='
    }
    energy_stats = {}
    for key, value in valid_pattern.iteritems():
        if value != '=':
            energy_stats[value] = None
    with open(det_filename) as f:
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
    with open(ct_filename) as f:
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


def process_fasta(in_fh, args, cluster_size_re, rna_seq_objs):
    """Process input file as fasta. Populate RNASequence (graph vertex)
    objects.
    """
    for record in SeqIO.parse(in_fh, 'fasta'):
        sequence = '%s%s%s'.replace('T', 'U') % (
            args.prefix, str(record.seq), args.suffix
        )
        cluster_size = 1
        try:
            cluster_size = cluster_size_re.search(record.description)
            cluster_size = cluster_size.group(1)
        except AttributeError:
            print 'Not able to find cluster size. Setting to 1.'
        if cluster_size is None:
            cluster_size = 1

        # find structure
        curr_seq = RNASequence(record.id, cluster_size, sequence)
        if args.run_mfold:
            curr_seq.structure, curr_seq.energy_dict = run_mfold(
                sequence, args
            )
            curr_seq.free_energy = curr_seq.energy_dict['dG']
        else:
            rnafold_out = run_rnafold(sequence, args)
            rnafold_out = rnafold_out.split('\n')
            try:
                curr_seq.structure, curr_seq.free_energy = (
                    rnafold_out[1].split(' (')
                )
            except (ValueError, IndexError):
                print 'Error running RNAfold:\n%s\nExiting.' % rnafold_out
                sys.exit(1)

            print '%s\n' % rnafold_out
            try:
                curr_seq.free_energy = abs(
                    float(curr_seq.free_energy.replace(')', ''))
                )
                curr_seq.ensemble_free_energy = abs(
                    float(rnafold_out[2].split('[')[1].replace(']', ''))
                )
                curr_seq.ensemble_probability = abs(float(
                    rnafold_out[4].split(';')[0].replace(
                        ' frequency of mfe structure in ensemble ', ''
                    )
                ))
                curr_seq.ensemble_diversity = abs(float(
                    rnafold_out[4].split(';')[1].replace(
                        ' ensemble diversity ', ''
                    )
                ))
            except IndexError:
                print (
                    'Error parsing RNAfold output. '
                    '(Couldn\'t find statistics.) Please check '
                    'RNAfold options.'
                )
                sys.exit(1)
        rna_seq_objs.append(curr_seq)


def process_struct_fasta(in_fh, args, cluster_size_re, rna_seq_objs):
    """Process non-fasta input file. Populate RNASequence
    (graph vertex) objects.
    """
    while True:
        try:
            # need to move through a triplet file structure, not fasta
            header, sequence, structure = list(
                itertools.islice(in_fh, 3)
            )
        except ValueError:  # end of file
            break
        sequence = sequence.strip('\n\r')
        structure = structure.strip('\n\r')
        if not structure.count('(') == structure.count(')'):
            continue
        if args.calc_structures:
            sequence = '%s%s%s'.replace('T', 'U') % (
                args.prefix, sequence, args.suffix
            )
        else:
            sequence = sequence.replace('T', 'U')
        try:
            cluster_size = cluster_size_re.search(header)
            cluster_size = cluster_size.group(1)
        except AttributeError:
            print 'Not able to find cluster size. Setting to 1.'
        if cluster_size is None:
            cluster_size = 1
        header = header.replace('>', '').strip()
        curr_seq = RNASequence(header, cluster_size, sequence)
        curr_seq.free_energy = 1
        curr_seq.ensemble_free_energy = 1
        curr_seq.ensemble_probability = 1
        curr_seq.ensemble_diversity = 1
        curr_seq.structure = structure
        rna_seq_objs.append(curr_seq)


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
        process_seq_pairs(seq_pairs, xgmml_obj, args, stats)
        for x in seq_pairs:
            if x.is_valid_edge:
                x.sequence2.use_for_comparison = False
        # delete nodes which already belong to an edge
        new_nodes_copy = [
            z for z in nodes_copy if z.use_for_comparison
        ]
        print 'Number of RNA sequences reduced from %d to %d ' % (
            len(nodes_copy), len(new_nodes_copy)
        )
        nodes_copy = new_nodes_copy

def find_edges_no_seed(rna_seq_objs, xgmml_obj, args, stats):
    """"Find edges using non-seed algorithm."""
    seq_pairs = []
    process_pairs_end = 0
    build_pairs_end = 0
    build_pairs_start = time.clock()

     # this makes the edges. looking at each pair of nodes
    for i in range(0, len(rna_seq_objs)):
        for j in range(i + 1, len(rna_seq_objs)):
            pair = RNASequencePair(
                rna_seq_objs[i], rna_seq_objs[j]
            )
            pair.energy_delta = abs(
                float(pair.sequence1.free_energy) - float(pair.sequence2.free_energy)
            )
            pair.edit_distance = Levenshtein.distance(
                pair.sequence1.sequence, pair.sequence2.sequence
            )
            seq_pairs.append(pair)
            # group things in batches of 10000  to find tree distances
            if len(seq_pairs) > 10000:
                build_pairs_end = time.clock()
                process_seq_pairs(seq_pairs, xgmml_obj, args, stats)
                procs_pairs_end = time.clock()

                build_time = build_pairs_end - build_pairs_start
                procs_time = procs_pairs_end - build_pairs_end

                print 'BUILD: {}'.format(build_time)
                print 'PROCS: {}'.format(procs_time)
                print 'PROCS > BUILD = {}'.format(procs_time > build_time)

                # zero out the seq_pairs array and start refilling again
                seq_pairs = []
                build_pairs_start = time.clock()

    # flush out the last of the tree distance seq_pairs
    process_seq_pairs(seq_pairs, xgmml_obj, args, stats)


def print_stats(stats, args):
    print '\nOverall Statistics:'
    print '--------------------'
    for stat in stats:
        stat_label = stat.capitalize().replace('_', ' ')
        if any(stats[stat]):
            print '%s mean: %.3g' % (stat_label, numpy.mean(stats[stat]))
            print '%s SD: %.3g' % (stat_label, numpy.std(stats[stat]))
            print '%s SEM: %.3g' % (stat_label, scipy.stats.sem(stats[stat]))
        else:
            print '%s mean: N/A' % (stat_label)
            print '%s SD: N/A' % (stat_label)
            print '%s SEM: N/A' % (stat_label)
        print

    for stat_pair in itertools.combinations(stats, 2):
        print 'Pearson correlation(%s, %s):' % (
            stat_pair[0], stat_pair[1]
        ),
        if args.calc_structures:
            pearson_coeff, p_value = scipy.stats.pearsonr(
                stats[stat_pair[0]], stats[stat_pair[1]]
            )
            print '%.3g, p=%.3g' % (pearson_coeff, p_value)
        else:
            print 'N/A'
