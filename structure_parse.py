#!/usr/bin/python
__author__ = 'tbair'


import sys
import re
import fileinput


split_re = re.compile('\)\.+\(')

#test_structure = "((((((((((...........))))))).((((((.((......))))))))...))).."

header = sys.stdin.readline()
sequence = sys.stdin.readline()
input_structure = sys.stdin.readline()
sub_structures = split_re.split(input_structure)
trimmed_structures = []

for structure in sub_structures:
    lefts = [m.start() for m in  re.finditer('\(',structure)]
    rights = [m.start() for m in re.finditer('\)',structure)]
    if len(lefts) >= len(rights):
        while len(lefts) >= len(rights):
            del lefts[0]
    elif (len(lefts) < len(rights)):
            while len(lefts) < len(rights)-1:
                    del rights[-1]

    trimmed_structures.append(structure[lefts[0]:rights[-1]])


trimmed_structures.sort(lambda x,y: cmp(len(x), len(y)))
for idx, trimmed_sub in enumerate(trimmed_structures):
    print header + str(idx+1)
    print trimmed_sub
