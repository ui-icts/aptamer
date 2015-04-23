#!/usr/bin/python
from Bio import SeqIO
import sys
import numpy
import subprocess
from Levenshtein import distance
import scipy.stats
from optparse import OptionParser
usage = "usage: %prog [options] /path/to/input/fasta > output_file"
parser = OptionParser(usage=usage)
parser.add_option("-v", dest="version", default=1, type="int", help= "Specifiy Vienna packager version 1 or 2 (1 is default)")
parser.add_option("-e", dest="editDistance", default=3, type="int", help= "Specifiy the minimum edit distance")
parser.add_option("-t", dest="treeDistance", default=3, type="int", help= "Specifiy the minimum tree distance")
parser.add_option("-n", dest="edgeType",  type="str", help= "Specify the type of edges to record (edit, tree, both)", default='both')
parser.add_option("-w", dest="write",  type="str", help= "Should or should not write xgml file (y or n) defaults to y", default='y')

(options,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit()
prefix = 'GGGAGGACGAUGCGG'
suffix = 'CAGACGACUCGCCCGA'
stats_energyDelta = []
stats_editDistance = []
stats_treeDistance = []
class Comparison:
    def __init__(self,S1,S2,xgmml):
        self.sequence1 = S1
        self.sequence2 = S2
        self.xgmml = xgmml
        self.energyDelta = None
        self.editDistance = None
        self.treeDistance = None
    def output(self):
        #output = [self.sequence1.name,
        #         self.sequence2.name,
        #         str(self.energyDelta),
        #         str(self.editDistance),
        #         str(self.treeDistance),
        #         ]
        #print ",".join(output)
        stats_energyDelta.append(self.energyDelta)
        stats_editDistance.append(self.editDistance)
        stats_treeDistance.append(float(self.treeDistance))
        if self.xgmml.nodes.has_key(self.sequence1.name) == False:
            self.xgmml.nodes[self.sequence1.name] = self.sequence1
        if self.xgmml.nodes.has_key(self.sequence2.name) == False:
            self.xgmml.nodes[self.sequence2.name] = self.sequence2
        if options.edgeType == 'both':
            if int(self.treeDistance) <= options.treeDistance:
                self.xgmml.edges.append([self.sequence1.name,self.sequence2.name,self.treeDistance,'treeDistance'])
            if int(self.editDistance) <= options.editDistance:
                self.xgmml.edges.append([self.sequence1.name,self.sequence2.name,self.treeDistance,'editDistance'])
        elif options.edgeType == 'edit':
            if int(self.editDistance) <= options.editDistance:
                self.xgmml.edges.append([self.sequence1.name,self.sequence2.name,self.treeDistance,'editDistance'])
        elif options.edgeType == 'tree':
            if int(self.treeDistance) <= options.treeDistance:
                self.xgmml.edges.append([self.sequence1.name,self.sequence2.name,self.treeDistance,'treeDistance'])
        else:
            print "Error in options %s not supported or recognized" % (options.editType,)
            parser.print_help()
            sys.exit()
            


class Sequence:
    def __init__(self,name,size,sequence):
        self.name = name
        self.clusterSize = size.replace('SIZE=','')
        self.sequence = sequence
        self.structure = None
        self.freeEnergy = None
        self.ensembleFreeEngergy = None
        self.ensembleProbablility = None
        self.ensembleDiversity = None
class XGMML:
    def __init__ (self,name):
        self.name  =  name 
        self.nodes = {}
        self.edges = []
    def output(self):
        string =  """<?xml version="1.0"?>\n<graph directed="1" id="5" label="%s"\nxmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\nxmlns:ns1="http://www.w3.org/1999/xlink"\nxmlns:dc="http://purl.org/dc/elements/1.1/"\nxmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\nxmlns="http://www.cs.rpi.edu/XGMML"> """ % (self.name,)
        string += "\n"
        for n in self.nodes:
            string +=  """<node id="%s" label="%s" weight="%s" >\n""" % (n,self.nodes[n].name,self.nodes[n].clusterSize)
            string += """<att type="integer" name="size" value="%s" />\n""" % (self.nodes[n].clusterSize,)
            string += """<att name="structure" value="%s" />\n""" % (self.nodes[n].structure,)
            string += """<att name="sequence" value="%s" /> \n""" % (self.nodes[n].sequence,)
            string += """<att name="energy" value="%s" /> \n""" % (self.nodes[n].freeEnergy,)
            string += """<att name="ensembleFreeEnergy" value="%s" /> \n""" % (self.nodes[n].ensembleFreeEnergy,)
            string += """<att name="ensembleProbability" value="%s" /> \n""" % (self.nodes[n].ensembleProbability,)
            string += """<att name="ensembleDiversity" value="%s" /> \n""" % (self.nodes[n].ensembleDiversity,)
            string += "</node>\n"
        for e in self.edges:
            string +=  """<edge source="%s" target="%s" label="%s to %s" >\n<att label="interaction" name="%s" value="%s" type="string"/>\n</edge>\n""" % (e[0],e[1],e[0],e[1],e[3],e[2])
            string += "\n"
        string +=  "</graph>\n"
        return string
def RNAFold(sequence,version):
    cmd = None
    if version == 1:
        cmd = ['RNAfold -p -T 30 -noLP -noPS -noGU']
    elif version == 2:
        cmd = ['RNAfold -p -T 30 --noLP --noPS --noGU']
    sffproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.PIPE,close_fds=True,shell=True)
    stdout_value, stderr_value = sffproc.communicate(sequence)
    return stdout_value

def RNAdistance(structures):
    cmd = ['RNAdistance']
    sffproc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.PIPE,close_fds=True,shell=True)
    stdout_value, stderr_value = sffproc.communicate(structures)
    return stdout_value

def mymean(a):
	return	float(sum(a))/len(a)

def processComparisons(comparisons):
    seqToTree = []
    for c in comparisons:
        seqToTree.append(c.sequence1.structure)
        seqToTree.append(c.sequence2.structure)
    treeDistance = RNAdistance("\n".join(seqToTree))
    treeDistance = treeDistance.strip("\n") #take off last lr
    treeDistance = treeDistance.split("\n")
    #assert len(treeDistance) == len(comparisons)
    for j in range(len(comparisons)):
        comparisons[j].treeDistance = treeDistance[j].split(' ')[1]  
        comparisons[j].output()

fastaHandle = open(args[0],'r')
data = []
#print "i,j,DiffEnergy,EditDist,TreeDist"
for record in SeqIO.parse(fastaHandle,"fasta"):
    sequence = prefix + str(record.seq) +suffix
    sequence = sequence.replace('T','U')
    thisseq = None
    try:
       thisseq = Sequence(record.id,record.description.split('-')[1],sequence)
    except:
       thisseq = Sequence(record.id,record.description.split('\t')[1],sequence)
    tmp = RNAFold(sequence,options.version)
    tmp = tmp.split("\n")
    thisseq.structure, thisseq.freeEnergy = tmp[1].split(' (')
    thisseq.freeEnergy = abs(float(thisseq.freeEnergy.replace(')','')))
    thisseq.ensembleFreeEnergy = abs(float(tmp[2].split('[')[1].replace(']','')))
    # frequency of mfe structure in ensemble 0.248667; ensemble diversity 8.19 
    thisseq.ensembleProbability = abs(float(tmp[4].split(';')[0].replace(' frequency of mfe structure in ensemble ','')))
    thisseq.ensembleDiversity= abs(float(tmp[4].split(';')[1].replace(' ensemble diversity ','')))
    data.append(thisseq)
xgmml = XGMML(sys.argv[1])
comparisons = []
for x in range(0, len(data)):
   for y in range (x+1, len(data)):
    comp = Comparison(data[x],data[y],xgmml)
    comparisons.append(comp)
    comp.energyDelta = abs(comp.sequence1.freeEnergy - comp.sequence2.freeEnergy)
    comp.editDistance = distance(comp.sequence1.sequence,comp.sequence2.sequence) #Lev edit distance
    if len(comparisons) > 10000:
        processComparisons(comparisons)
        comparisons = []
processComparisons(comparisons)

if options.write == 'y':
    cytoscapeOut = open (args[0]+".xgmml",'w')
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
print "%s %s pearsons corr tree:edit" % scipy.stats.pearsonr(stats_treeDistance,stats_editDistance)
print "%s %s pearsons corr tree:energy" % scipy.stats.pearsonr(stats_treeDistance,stats_energyDelta)
print "%s %s pearsons corr edit:energy" % scipy.stats.pearsonr(stats_editDistance,stats_energyDelta)

#stats_energyDelta = []
#stats_editDistance = []
#stats_treeDistance = []

