#!/usr/bin/python
from Bio import SeqIO
import os
import sys
import subprocess
from Levenshtein import distance
import scipy.stats
import re
import numpy
from itertools import islice
from optparse import OptionParser


usage = "usage: %prog [options] /path/to/input/fasta > output_file"

parser = OptionParser(usage=usage)

parser.add_option("-a", dest="seed", action="store_true", default=False,
                  help="Specify this if you want to use the seed sequence algo")

parser.add_option("-e", dest="editDistance", default=3, type="int", help="Specify the minimum edit distance")

parser.add_option("-f", dest="fasta", action="store_true", default=False,
                  help="Set this if the input contains a structure line")

parser.add_option("-m", dest="runMfold", action="store_true", default=False,
                  help="Specify whether to add Mfold structures.")

parser.add_option("-n", dest="edgeType", type="str", help="Specify the type of edges to record (edit, tree, both)",
                  default='both')

parser.add_option("-p", dest="prefix", type="str",
                  help="The prefix sequence used for structure prediction defaults to GGGAGGACGAUGCG",
                  default='GGGAGGACGAUGCGG')

parser.add_option("-r", dest="rna", action="store_true", default=False,
                  help="Specify this if your sequence is RNA ie has U instead of T")

parser.add_option("-s", dest="suffix", type="str",
                  help="The prefix sequence used for structure prediction defaults to CAGACGACUCGCCCGA",
                  default='CAGACGACUCGCCCGA')

parser.add_option("-t", dest="treeDistance", default=3, type="int", help="Specify the minimum tree distance")

parser.add_option("-v", dest="version", default=2, type="int",
                  help="Specify Vienna packager version 1 or 2 (1 is default)")

parser.add_option("-w", dest="write", type="str", help="Should or should not write xgmml file (y or n) defaults to y",
                  default='y')

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.print_help()
    sys.exit()
sizere = re.compile("SIZE=(\d+)")  # set up re for the cluster size
fastaHandle = open(args[0], 'r')
# use this if the input is in fasta format
stats_energyDelta = []
stats_editDistance = []
stats_treeDistance = []


class Comparison:
    def __init__(self, s1, s2, xgmml_pointer):
        self.sequence1 = s1
        self.sequence2 = s2
        self.xgmml = xgmml_pointer
        self.energyDelta = None
        self.editDistance = None
        self.treeDistance = None
        self.flag = False

    def matched(self):
        if self.flag:
            return True
        else:
            return False

    def output(self):
        stats_energyDelta.append(self.energyDelta)
        stats_editDistance.append(self.editDistance)
        stats_treeDistance.append(float(self.treeDistance))
        if not self.sequence1.name in self.xgmml.nodes:  # if the xgmml data structure does not have this node add it
            self.xgmml.nodes[self.sequence1.name] = self.sequence1
        if not self.sequence2.name in self.xgmml.nodes:  # if the xgmml data structure does not have this node add it
            self.xgmml.nodes[self.sequence2.name] = self.sequence2

        if options.edgeType == 'both':
            if int(self.treeDistance) <= options.treeDistance:
                self.xgmml.edges.append([self.sequence1.name, self.sequence2.name, self.treeDistance, 'treeDistance'])
                self.flag = True
            if int(self.editDistance) <= options.editDistance:
                self.xgmml.edges.append([self.sequence1.name, self.sequence2.name, self.editDistance, 'editDistance'])
                self.flag = True  # have a match need to remove these from future consideration
        elif options.edgeType == 'edit':
            if int(self.editDistance) <= options.editDistance:
                self.xgmml.edges.append([self.sequence1.name, self.sequence2.name, self.treeDistance, 'editDistance'])
                self.flag = True  # have a match need to remove these from future consideration
        elif options.edgeType == 'tree':
            if int(self.treeDistance) <= options.treeDistance:
                self.xgmml.edges.append([self.sequence1.name, self.sequence2.name, self.treeDistance, 'treeDistance'])
                self.flag = True  # have a match need to remove these from future consideration
        else:
            print "Error in options %s not supported or recognized" % (options.editType,)
            parser.print_help()
            sys.exit()


class Sequence:
    def __init__(self, name, si, seq):
        self.name = name
        self.clusterSize = int(si)
        self.sequence = seq
        self.structure = None
        self.freeEnergy = None
        self.ensembleFreeEnergy = None
        self.ensembleProbability = None
        self.ensembleDiversity = None
        self.useForComparison = True

    def full_output(self):
        attrs = vars(self)
        print ','.join("%s:%s" % item for item in attrs.items())

    def output(self):
        print ">%s  SIZE=%s" % (self.name, self.clusterSize)
        print "%s" % (self.sequence,)
        print "%s" % (self.structure,)


class XGMML:
    def __init__(self, name):
        self.name = name
        self.nodes = {}
        self.edges = []

    def output(self):
        string = """<?xml version="1.0"?>\n<graph directed="1" id="5" label="%s"\n\
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n\
        xmlns:ns1="http://www.w3.org/1999/xlink"\nxmlns:dc="http://purl.org/dc/elements/1.1/"\n\
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\nxmlns="http://www.cs.rpi.edu/XGMML"> """ % (
            self.name,)
        string += "\n"
        for n in self.nodes:
            string += """<node id="%s" label="%s" weight="%s" >\n""" % (
                n, self.nodes[n].name, self.nodes[n].clusterSize)
            string += """<att type="integer" name="size" label="Size" value="%s" />\n""" % (self.nodes[n].clusterSize,)
            string += """<att type="string" label="Structure" name="structure" value="%s" />\n""" % (self.nodes[n].structure,)
            if options.runMfold:
                string += """<att type="string" label="mfold Structure" name="mfoldStructure" value="%s" />\n""" % (self.nodes[n].mfoldStructure,)
            string += """<att type="string" label="Sequence" name="sequence" value="%s" /> \n""" % (self.nodes[n].sequence,)
            string += """<att type="real" label="Energy" name="energy" value="%s" /> \n""" % (self.nodes[n].freeEnergy,)
            string += """<att type="real" label="ensemble Free Energy" name="ensembleFreeEnergy" value="%s" /> \n""" % (self.nodes[n].ensembleFreeEnergy,)
            string += """<att type="real" label="ensemble Probability" name="ensembleProbability" value="%s" /> \n""" % (self.nodes[n].ensembleProbability,)
            string += """<att type="real" label="ensemble Diversity" name="ensembleDiversity" value="%s" /> \n""" % (self.nodes[n].ensembleDiversity,)
            string += "</node>\n"
        for e in self.edges:
            string += """<edge source="%s" target="%s" label="%s to %s" >\n\
            <att label="interaction" name="%s" value="%s" type="string"/>\n</edge>\n""" % (
                e[0], e[1], e[0], e[1], e[3], e[2])
            string += "\n"
        string += "</graph>\n"
        return string


def main():
    data = []  # main data structure is a list of Sequence objects
    if options.fasta:
        for record in SeqIO.parse(fastaHandle, "fasta"):
            sequence = options.prefix + str(record.seq) + options.suffix
            if not options.rna:
                sequence = sequence.replace('T', 'U')
            size = 1

            try:
                size = sizere.search(record.description)
                size = size.group(1)
            except AttributeError:
                print 'Not able to find size, setting to 1'

            thisseq = Sequence(record.id, size, sequence)
            tmp = RNAFold(sequence, options.version)
            tmp = tmp.split("\n")
            print tmp
            thisseq.structure, thisseq.freeEnergy = tmp[1].split(' (')
            thisseq.freeEnergy = abs(float(thisseq.freeEnergy.replace(')', '')))
            thisseq.ensembleFreeEnergy = abs(float(tmp[2].split('[')[1].replace(']', '')))
            # frequency of mfe structure in ensemble 0.248667; ensemble diversity 8.19
            thisseq.ensembleProbability = abs(
                float(tmp[4].split(';')[0].replace(' frequency of mfe structure in ensemble ', '')))
            thisseq.ensembleDiversity = abs(float(tmp[4].split(';')[1].replace(' ensemble diversity ', '')))
            data.append(thisseq)
            if options.runMfold:
                mfold_structure = run_mfold(sequence)
                thisseq.mfoldStructure = mfold_structure
            else:
                thisseq.mfoldStructure = None

    else:  # use this if the data is not in fasta format
        while True:
            try:
                header, sequence, structure = list(
                    islice(fastaHandle, 3))  # need to move through a triplet file structure not fasta
            except ValueError:
                break
            sequence = sequence.strip('\n')
            sequence = sequence.strip('\r')
            structure = structure.strip('\n')
            structure = structure.strip('\r')
            if not structure.count('(') == structure.count(')'):
                continue
            sequence = options.prefix + sequence + options.suffix
            if not options.rna:
                sequence = sequence.replace('T', 'U')
            size = 1
            try:
                size = sizere.search(header)
                size = size.group(1)
            except AttributeError:
                print "Not able to find size setting to 1"
            header = header.replace('>', '')
            header = header.split('SIZE=')[0]
            thisseq = Sequence(header, size, sequence)
            thisseq.freeEnergy = 1
            thisseq.ensembleFreeEnergy = 1
            thisseq.ensembleProbability = 1
            thisseq.ensembleDiversity = 1
            thisseq.structure = structure
            data.append(thisseq)

    xgmml = XGMML(sys.argv[1])
    #data structure should now be populated, now need to move through and find all the connections
    if options.seed:
        past_harvested_nodes = {}  # only used for seed algo but collected in all cases
        while len(data) > 2:
            data[0].useForComparison = False
            comparisons = []
             # go through and find all the matches, then remove matches and start again till gone
            for x in range(1, len(data)-1):
                comp = Comparison(data[0], data[x], xgmml)
                comp.energyDelta = abs(comp.sequence1.freeEnergy - comp.sequence2.freeEnergy)
                # distance function imported from from Levenshtein
                comp.editDistance = distance(comp.sequence1.sequence, comp.sequence2.sequence)
                comparisons.append(comp)
            processComparisons(comparisons)
            for c in comparisons:
                if c.matched():
                    c.sequence2.useForComparison = False
            newData = []
            for d in data:
                if d.useForComparison:
                    newData.append(d)
            print "%s reduced to %s " % (len(data), len(newData))
            data = newData

    else:
        comparisons = []
        for x in range(0, len(data)):  # this makes the edges
            for y in range(x + 1, len(data)):
                comp = Comparison(data[x], data[y], xgmml)
                comp.energyDelta = abs(comp.sequence1.freeEnergy - comp.sequence2.freeEnergy)
                # distance function imported from from Levenshtein
                comp.editDistance = distance(comp.sequence1.sequence, comp.sequence2.sequence)
                comparisons.append(comp)
                if len(comparisons) > 10000:  # group things in batches of 10000  to find tree distances
                    processComparisons(comparisons)
                    comparisons = []  # zero out the comparisons array and start refilling again
        processComparisons(comparisons)  # flush out the last of the tree distance comparisons

    if options.write == 'y':
        cytoscapeOut = open(args[0] + ".xgmml", 'w')
        cytoscapeOut.write(xgmml.output())
        cytoscapeOut.close()
    print "%s mean energyDelta" % numpy.mean(stats_energyDelta)
    print "%s std energyDelta" % numpy.std(stats_energyDelta)
    print "%s sem energyDelta" % scipy.stats.sem(stats_energyDelta)
    print "%s mean editDistance" % numpy.mean(stats_editDistance)
    print "%s std editDistance" % numpy.std(stats_editDistance)
    print "%s sem editDistance" % scipy.stats.sem(stats_editDistance)
    print "%s mean treeDistance" % numpy.mean(stats_treeDistance)
    print "%s std treeDistance" % numpy.std(stats_treeDistance)
    print "%s sem treeDistance" % scipy.stats.sem(stats_treeDistance)
    if options.fasta:
        print "%s %s pearsons corr tree:edit" % scipy.stats.pearsonr(stats_treeDistance, stats_editDistance)
        print "%s %s pearsons corr tree:energy" % scipy.stats.pearsonr(stats_treeDistance, stats_energyDelta)
        print "%s %s pearsons corr edit:energy" % scipy.stats.pearsonr(stats_editDistance, stats_energyDelta)
    else:  # ugly hack because we don't have the energy values for things we didn't calculate
        print "1 1 pearsons corr tree:edit"
        print "1 1 pearsons corr tree:energy"
        print "1 1 pearsons corr edit:energy"


def RNAFold(seq, version):
    print 'Running RNAFold...'
    cmd = None
    if version == 1:
        cmd = ['RNAfold -p -T 30 -noLP -noPS -noGU']
    elif version == 2:
        cmd = ['RNAfold -p -T 30 --noLP --noPS --noGU']
    sffproc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.PIPE,
                               close_fds=True, shell=True)
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
    x = subprocess.Popen(
        ['mfold SEQ=%s T=30' % temp_filename], stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, close_fds=True, shell=True
    )
    stdout, stderr = x.communicate()
    mfold_structure = convert_ct_to_bracket_dot('%s.ct' % temp_filename)
    os.chdir('..')
    return mfold_structure


def RNAdistance(structures):
    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.PIPE,
                               close_fds=True, shell=True)
    stdout_value, stderr_value = sffproc.communicate(structures)
    return stdout_value


def mymean(a):
    return float(sum(a)) / len(a)


def processComparisons(comparison_list):
    #helper function that runs a batch of comparisons
    seqToTree = []
    for c in comparison_list:
        seqToTree.append(c.sequence1.structure)
        seqToTree.append(c.sequence2.structure)
    treeDistance = RNAdistance("\n".join(seqToTree))
    treeDistance = treeDistance.strip("\n")  # take off last lr
    treeDistance = treeDistance.split("\n")
    assert len(treeDistance) == len(comparison_list)
    for j in range(len(comparison_list)):
        comparison_list[j].treeDistance = treeDistance[j].split(' ')[1]
        comparison_list[j].output()


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
    main()
