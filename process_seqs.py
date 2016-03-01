#!/usr/bin/python
import os
import sys
import re
import subprocess
import itertools
import textwrap
import optparse
import numpy
import scipy.stats
from Bio import SeqIO
import Levenshtein


class RNASequence:
    def __init__(self, name, cluster_size, seq):
        self.name = name
        self.cluster_size = int(cluster_size)
        self.sequence = seq
        self.structure = None
        self.free_energy = None
        self.ensemble_free_energy = None
        self.ensemble_probability = None
        self.ensemble_diversity = None
        self.use_for_comparison = True

    def full_output(self):
        attrs = vars(self)
        print ','.join('%s:%s' % item for item in attrs.items())

    def output(self):
        print '>%s  SIZE=%s' % (self.name, self.cluster_size)
        print self.sequence
        print self.structure


class RNASequencePair:
    """Pair of RNASequence objects"""
    def __init__(self, seq1, seq2, xgmml_obj):
        self.sequence1 = seq1
        self.sequence2 = seq2
        self.xgmml = xgmml_obj
        self.energy_delta = None
        self.edit_distance = None
        self.tree_distance = None
        self.match_exists = False

    def output(self, options):
        # if the xgmml data structure does not have this node, add it
        if not self.sequence1.name in self.xgmml.nodes:  
            self.xgmml.nodes[self.sequence1.name] = self.sequence1
        if not self.sequence2.name in self.xgmml.nodes:
            self.xgmml.nodes[self.sequence2.name] = self.sequence2

        if options.edge_type == 'both':
            if (
                self.tree_distance and
                (int(self.tree_distance) <= options.tree_distance)
            ):
                self.xgmml.edges.append([
                    self.sequence1.name, self.sequence2.name,
                    self.tree_distance, 'tree_distance'
                ])
                self.match_exists = True
            if (
                self.edit_distance and
                (int(self.edit_distance) <= options.edit_distance)
            ):
                self.xgmml.edges.append([
                    self.sequence1.name, self.sequence2.name,
                    self.edit_distance, 'edit_distance'
                ])
                # have a match. need to remove these from future consideration
                self.match_exists = True
        elif options.edge_type == 'edit':
            if (
                self.edit_distance and
                (int(self.edit_distance) <= options.edit_distance)
            ):
                self.xgmml.edges.append([
                    self.sequence1.name, self.sequence2.name,
                    self.tree_distance, 'edit_distance'
                ])
                self.match_exists = True
        elif options.edge_type == 'tree':
            if (
                self.tree_distance and
                (int(self.tree_distance) <= options.tree_distance)
            ):
                self.xgmml.edges.append([
                    self.sequence1.name, self.sequence2.name,
                    self.tree_distance, 'tree_distance'
                ])
                self.match_exists = True
        else:
            print 'Error: Option not supported or recognized: %s' % (
                options.edge_type
            )
            parser.print_help()
            sys.exit()


class XGMML:
    """XGMML file"""
    def __init__(self, name):
        self.name = name
        self.nodes = {}
        self.edges = []
        self.out_str = ''

    def output_att(self, type_, name, label, value):
        self.out_str += '<att type="%s" name="%s" label="%s" value="%s">' % (
            type_, name, label, value
        )

    def output(self, options):
        self.out_str = textwrap.dedent(
            """
            <?xml version="1.0"?>
            <graph directed="1" id="5" label="%s"
            xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
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
            self.output_att(
                'integer', 'size', 'Size', self.nodes[n].cluster_size
            )
            self.output_att(
                'string', 'structure', 'Structure', self.nodes[n].structure
            )
            if options.run_mfold:
                self.output_att(
                    'string', 'mfoldStructure', 'mfold Structure',
                    self.nodes[n].mfoldStructure
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
            self.out_str += (
                '<edge source="%s" target="%s" label="%s to %s" >\n' % (
                    e[0], e[1], e[0], e[1]
                )
            )
            self.output_att('string', e[3], 'interaction', e[2])
            self.out_str += '</edge>\n\n'
        self.out_str += '</graph>\n'
        return self.out_str


####### MAIN #######
def main():
    options, args = parse_arguments()

    cluster_size_re = re.compile('SIZE=(\d+)')
    fasta_fh = open(args[0])

    # stats to be used if the input is in fasta format (fasta option selected)
    # stats are given default values if input is non-fasta
    stats = {'energy_delta':[], 'edit_distance':[], 'tree_distance':[]}

    data = []  # list of Sequence objects
    if options.fasta:
        for record in SeqIO.parse(fasta_fh, 'fasta'):
            sequence = '%s%s%s'.replace('T', 'U') % (
                options.prefix, str(record.seq), options.suffix
            )
            cluster_size = 1

            try:
                cluster_size = cluster_size_re.search(record.description)
                cluster_size = cluster_size.group(1)
            except AttributeError:
                print 'Not able to find cluster size. Setting to 1.'

            curr_seq = RNASequence(record.id, cluster_size, sequence)
            rnafold_out = run_rnafold(sequence, options.vienna_version)
            rnafold_out = rnafold_out.split('\n')
            print rnafold_out
            curr_seq.structure, curr_seq.free_energy = (
                rnafold_out[1].split(' (')
            )
            curr_seq.free_energy = abs(
                float(curr_seq.free_energy.replace(')', ''))
            )
            curr_seq.ensemble_free_energy = abs(
                float(rnafold_out[2].split('[')[1].replace(']', ''))
            )
            # frequency of mfe structure in ensemble 0.248667;
            # ensemble diversity 8.19
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
            data.append(curr_seq)
            if options.run_mfold:
                mfold_structure = run_mfold(sequence)
                curr_seq.mfoldStructure = mfold_structure
            else:
                curr_seq.mfoldStructure = None

    else:  # use this if the data is not in fasta format
        while True:
            try:
                # need to move through a triplet file structure, not fasta
                header, sequence, structure = list(
                    itertools.islice(fasta_fh, 3)
                )
            except ValueError:
                break
            sequence = sequence.strip('\n')
            sequence = sequence.strip('\r')
            structure = structure.strip('\n')
            structure = structure.strip('\r')
            if not structure.count('(') == structure.count(')'):
                continue
            sequence = '%s%s%s'.replace('T', 'U') % (
                options.prefix, sequence, options.suffix
            )
            try:
                cluster_size = cluster_size_re.search(header)
                cluster_size = cluster_size.group(1)
            except AttributeError:
                print 'Not able to find cluster size. Setting to 1.'
            if cluster_size is None:
                cluster_size = 1
            header = header.replace('>', '')
            header = header.split('SIZE=')[0]
            curr_seq = RNASequence(header, cluster_size, sequence)
            curr_seq.free_energy = 1
            curr_seq.ensemble_free_energy = 1
            curr_seq.ensemble_probability = 1
            curr_seq.ensemble_diversity = 1
            curr_seq.structure = structure
            data.append(curr_seq)

    xgmml = XGMML(sys.argv[1])
    #data structure should now be populated, now need to move through and find all the connections
    if options.seed:
        past_harvested_nodes = {}  # only used for seed algo but collected in all cases
        while len(data) > 2:
            data[0].use_for_comparison = False
            comparisons = []
             # go through and find all the matches, then remove matches and start again till gone
            for x in range(1, len(data)-1):
                comp = RNASequencePair(data[0], data[x], xgmml)
                comp.energy_delta = abs(
                    comp.sequence1.free_energy - comp.sequence2.free_energy
                )
                comp.edit_distance = Levenshtein.distance(
                    comp.sequence1.sequence, comp.sequence2.sequence
                )
                comparisons.append(comp)
            process_comparisons(comparisons, options, stats)
            for c in comparisons:
                if c.match_exists:
                    c.sequence2.use_for_comparison = False
            newData = []
            for d in data:
                if d.use_for_comparison:
                    newData.append(d)
            print '%s reduced to %s ' % (len(data), len(newData))
            data = newData

    else:
        comparisons = []
        for x in range(0, len(data)):  # this makes the edges
            for y in range(x + 1, len(data)):
                comp = RNASequencePair(data[x], data[y], xgmml)
                comp.energy_delta = abs(
                    comp.sequence1.free_energy - comp.sequence2.free_energy
                )
                comp.edit_distance = Levenshtein.distance(
                    comp.sequence1.sequence, comp.sequence2.sequence
                )
                comparisons.append(comp)
                if len(comparisons) > 10000:  # group things in batches of 10000  to find tree distances
                    process_comparisons(comparisons, options, stats)
                    comparisons = []  # zero out the comparisons array and start refilling again
        process_comparisons(comparisons, options, stats)  # flush out the last of the tree distance comparisons

    if not options.no_xgmml:
        cytoscapeOut = open(args[0] + '.xgmml', 'w')
        cytoscapeOut.write(xgmml.output(options))
        cytoscapeOut.close()

    print '\n\nOverall Statistics:'
    for stat in stats:
        stat_label = stat.capitalize().replace('_', ' ')
        if any(stats[stat]):
            print '%s mean: %.2f' % (stat_label, numpy.mean(stats[stat]))
            print '%s SD: %.2f' % (stat_label, numpy.std(stats[stat]))
            print '%s SEM: %.2f' % (stat_label, scipy.stats.sem(stats[stat]))
        else:
            print '%s mean: N/A' % (stat_label)
            print '%s SD: N/A' % (stat_label)
            print '%s SEM: N/A' % (stat_label)
        print

    for stat_pair in itertools.combinations(stats, 2):
        print 'Pearson correlation(%s, %s):' % (
            stat_pair[0], stat_pair[1]
        ),
        if options.fasta:
            pearson_coeff, p_value = scipy.stats.pearsonr(
                stats[stat_pair[0]], stats[stat_pair[1]]
            )
            print '%.2f, p=%.3f' % (pearson_coeff, p_value)
        else:
            print 'N/A'

    fasta_fh.close()


def parse_arguments():
    parser = optparse.OptionParser(
        usage='usage: %prog [options] /path/to/input/fasta > output_file'
    )

    parser.add_option(
        '-a', dest='seed', action='store_true', default=False,
        help='Specify this if you want to use the seed sequence algo'
    )
    parser.add_option(
        '-e', dest='edit_distance', default=3, type='int',
        help='Specify the minimum edit distance'
    )
    ##-f just means that there's a SIZE= in input file?
    parser.add_option(
        '-f', dest='fasta', action='store_true', default=False,
        help='Set this if the input contains a structure line'
    )
    parser.add_option(
        '-m', dest='run_mfold', action='store_true', default=False,
        help='Specify whether to add Mfold structures.'
    )
    parser.add_option(
        '-n', dest='edge_type', type='str',
        help='Specify the type of edges to record (edit, tree, both)',
        default='both'
    )
    parser.add_option(
        '-p', dest='prefix', type='str',
         help=(
            'The prefix sequence used for structure prediction. '
            'Default is GGGAGGACGAUGCG'
        ),
        default='GGGAGGACGAUGCGG'
    )
    parser.add_option(
        '-s', dest='suffix', type='str',
        help= (
            'The prefix sequence used for structure prediction. '
            'Default is CAGACGACUCGCCCGA'
        ),
        default='CAGACGACUCGCCCGA'
    )
    parser.add_option(
        '-t', dest='tree_distance', default=3, type='int',
        help='Specify the minimum tree distance'
    )
    parser.add_option(
        '-v', dest='vienna_version', default=2, type='int',
        help='Specify Vienna packager version 1 or 2 (1 is default)'
    )
    parser.add_option(
        '-x', dest='no_xgmml', action='store_true', default=False,
        help='XGMML output file will NOT be created if this '
        'option is specified'
    )

    options, args = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    return (options, args)


def run_rnafold(seq, vienna_version):
    print 'Running RNAFold...'
    cmd = None
    if vienna_version == 1:
        cmd = ['RNAfold -p -T 30 -noLP -noPS -noGU']
    elif vienna_version == 2:
        cmd = ['RNAfold -p -T 30 --noLP --noPS --noGU']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True
    )
    stdout_value, stderr_value = sffproc.communicate(seq)
    return stdout_value


def run_mfold(seq):
    print 'Running mfold...'
    if not os.path.exists('mfold_out'):
        os.mkdir('mfold_out')
    os.chdir('mfold_out')
    temp_filename = 'mfold_temp.txt'
    with open(temp_filename, 'w') as f:
        f.write('%s\n' % seq)

    ret = subprocess.Popen(
        ['mfold SEQ=%s T=30' % temp_filename], stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, close_fds=True, shell=True
    )
    stdout, stderr = ret.communicate()
    if ret.returncode != 0:
        mfold_error = ''
        if stdout:
            mfold_error += str(stdout)
        if stderr:
            mfold_error += str(stderr)
        print 'Error when running mfold:\n%s' % (mfold_error)
        sys.exit(ret.returncode)
    mfold_structure = convert_ct_to_bracket_dot('%s.ct' % temp_filename)
    os.chdir('..')
    return mfold_structure


def rna_distance(structures):
    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True
    )
    stdout_value, stderr_value = sffproc.communicate(structures)
    return stdout_value


def process_comparisons(comparison_list, options, stats):
    #helper function that runs a batch of comparisons
    seqToTree = []
    for c in comparison_list:
        seqToTree.append(c.sequence1.structure)
        seqToTree.append(c.sequence2.structure)
    tree_distance = rna_distance('\n'.join(seqToTree))
    tree_distance = tree_distance.strip('\n')  # take off last lr
    tree_distance = tree_distance.split('\n')
    assert len(tree_distance) == len(comparison_list)
    for i, c in enumerate(comparison_list):
        comparison_list[i].tree_distance = tree_distance[i].split(' ')[1]
        c.output(options)
        stats['energy_delta'].append(c.energy_delta)
        stats['edit_distance'].append(c.edit_distance)
        try:
            stats['tree_distance'].append(float(c.tree_distance))
        except ValueError:
            stats['tree_distance'].append(None)


def convert_ct_to_bracket_dot(ct_filename):
    bracket_dot = ''
    with open(ct_filename) as f:
        for row in f:
            row = row.split()
            if '=' in row and len(row) < 6:  # first row
                energy = row[3]
            elif row[4] == '0':
                bracket_dot += '.'
            elif int(row[0]) < int(row[4]):
                bracket_dot += '('
            elif int(row[0]) > int(row[4]):
                bracket_dot += ')'
    return '%s (%s)' % (bracket_dot, energy) if (bracket_dot != '') else None


if __name__ == '__main__':
    sys.exit(main())
