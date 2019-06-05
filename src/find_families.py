"""Find connected components from an xgmml file"""
from __future__ import print_function

from builtins import str
from builtins import range
from builtins import object
import sys
import os
import argparse

class Edge(object):
    """Graph edge."""
    def __init__(self, source=None, target=None, edit_dist=None):
        self.source = source  # name of source node
        self.target = target  # name of target node
        self.edit_dist = edit_dist

    def __repr__(self):
        return str(self.__dict__)


def main():
    args = parse_arguments()
    with open(args.input_file) as in_f:
        nodes, edges = read_in_input_file(in_f)

    families = {}  # {max edit dist: list of connected components}
    for max_edit_dist in range(1, 5):  # max edit distance 1, 2, 3, and >=4
        # {node: list of connected nodes}
        conn_nodes = find_connected_nodes(nodes, edges, max_edit_dist)

        # do depth-first search on each unvisited node
        # to find other nodes in family
        visited_nodes = []
        families[max_edit_dist] = []
        for node in conn_nodes:
            if node in visited_nodes:
                continue
            family = [node]
            visited_nodes.append(node)
            dfs(node, conn_nodes, visited_nodes, family)  # updates family
            families[max_edit_dist].append(family)

    # output
    if args.output_file:
        out_f_name = args.output_file
    else:
        out_f_name = '%s.tsv' % args.input_file
    with open(out_f_name, 'w') as out_f:
        output_families(nodes, families, out_f)
    print('Output written to %s' % out_f_name)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            '%s: Output families (connected components) for ' 
            'each edit distance.' % os.path.basename(sys.argv[0])
        )
    )
    parser.add_argument(
        'input_file', help='Specify path of input xgmml graph file.'
    )
    parser.add_argument(
        '-o', '--output_file', help=(
            'Specify path of output tab-separated file. '
            '(Default: <input_filename>.tsv)'
        )
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        print('\nError: Input file not specified.')
        sys.exit(1)

    return parser.parse_args()


def find_connected_nodes(nodes, edges, max_edit_dist):
    """Return dict of node : list of connected nodes"""
    conn_nodes = {}
    for edge in edges:
        if (edge.edit_dist <= max_edit_dist or max_edit_dist == 4):
            if edge.source not in conn_nodes:
                conn_nodes[edge.source] = []
            conn_nodes[edge.source].append(edge.target)
            if edge.target not in conn_nodes:
                conn_nodes[edge.target] = []
            conn_nodes[edge.target].append(edge.source)
    for node in nodes:
        if node not in conn_nodes:
            conn_nodes[node] = []
    return conn_nodes


def dfs(node, conn_nodes, visited_nodes, family):
    for n in conn_nodes[node]:
        if n in visited_nodes:
            continue
        family.append(n)
        visited_nodes.append(n)
        dfs(n, conn_nodes, visited_nodes, family)


def read_in_input_file(in_f):
    """Read in and validate nodes and edges from input xgmml."""
    nodes = []
    edges = []
    source = None
    target = None
    edit_dist = None

    for line in in_f:
        line = line.strip()
        if line.startswith('<node'):
            nodes.append(line.split('"')[1].split('"')[0])
        elif line.startswith('<edge'):
            source = line.split('source="')[1].split('"')[0]
            target = line.split('target="')[1].split('"')[0]
        # looking for specific string 'editDistance'
        elif 'editDistance' in line:
            edit_dist = int(line.split('value="')[1].split('"')[0])

        elif line.startswith('</edge'):
            if not source or not target:
                print('Could not find source/ target for edge.')
                sys.exit(1)
            if edit_dist is None:
                if source and target:
                    print (
                        'Could not find edit distance for '
                        'edge (%s, %s).' % (source, target)
                    )
                else:
                    print('Could not find edit distance for edge.')
                sys.exit(1)
            edges.append(
                Edge(source=source, target=target, edit_dist=edit_dist)
            )
            source = None
            target = None
            edit_dist = None

    # validate nodes, edges
    if not nodes:
        print (
            'Could not find any nodes in input file. Please check '
            'that file is in xgmml format.'
        )
        sys.exit(1)
    if not edges:
        print (
            'Could not find any edges in input file. Please check '
            'that file is in xgmml format.'
        )
        sys.exit(1)

    return (nodes, edges)


def output_families(nodes, families, out_f):
    edit_dist_strs = ['=1', '=2', '=3', '>=4']
    out_f.write(
        '\t%s\n' % '\t'.join(['MaxEditDist%s' % z for z in edit_dist_strs])
    )
    # {(node, edit distance): family}}
    d = {}
    for edit_dist in families:
        # i is family name
        for i, family in enumerate(families[edit_dist]):
            for node in family:
                d[(node, edit_dist)] = i
    #try:
    #    nodes = [int(z) for z in nodes]
    #except ValueError:
    #    pass
    for node in sorted(nodes):
        node = str(node)
        out_f.write('Node=%s\t' % node)
        out_f.write('%s\n' % '\t'.join([
            'Family=%s' % d[(node, edit_dist)] for edit_dist in sorted(families)
        ]))


if __name__ == main():
    sys.exit(main())
