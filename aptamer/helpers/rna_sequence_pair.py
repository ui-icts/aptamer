import RNA
import Levenshtein
from helpers import functions


class RNASequencePair:

    __slots__ = 'sequence1', 'sequence2', 'tree_distance', 'energy_delta', 'edit_distance', 'is_valid_edge'

    @staticmethod
    def build(pair):
        seq1, seq2 = pair
        td = RNASequencePair.compute_tree_distance(seq1,seq2)
        return RNASequencePair(seq1,seq2,td)

    @staticmethod
    def compute_tree_distance(sequence1,sequence2):
        t1 = RNA.make_tree(RNA.expand_Full(sequence1.structure))
        t2 = RNA.make_tree(RNA.expand_Full(sequence2.structure))
        value = RNA.tree_edit_distance(t1, t2)
        RNA.free_tree(t1)
        RNA.free_tree(t2)
        return str(int(value))

    """Graph edge. Pair of RNASequence objects."""
    def __init__(self, seq1, seq2, tree_distance=None):
        self.sequence1 = seq1
        self.sequence2 = seq2
        self.tree_distance = tree_distance

        self.energy_delta = None
        self.edit_distance = None
        self.is_valid_edge = False  # used only with seed option

        self.update_energy_delta()
        self.update_edit_distance()

    def __str__(self):
        return '%s\n---\n%s' % (str(self.sequence1), str(self.sequence2))

    def update_energy_delta(self):
        fe1 = float(self.sequence1.free_energy)
        fe2 = float(self.sequence2.free_energy)
        self.energy_delta = abs(fe1 - fe2)

    def update_edit_distance(self):
        self.edit_distance = Levenshtein.distance(
                self.sequence1.sequence,
                self.sequence2.sequence
                                                )

    def update_tree_distance(self):

        seq_to_tree = []
        seq_to_tree.append(self.sequence1.structure)
        seq_to_tree.append(self.sequence2.structure)
        string_of_sequences = '\n'.join(seq_to_tree)

        tree_distance = functions.rna_distance(string_of_sequences)
        # take off last lr
        tree_distance = tree_distance.strip('\n').split('\n')
        assert len(tree_distance) == 1, (
            'Error length of tree distance %s does not match length of seq_pairs '
            '%s and should -- check installation of RNAdistance'
            % (len(tree_distance), 1)
        )

        self.tree_distance = tree_distance[0].split(' ')[1]
        # test_value = self.compute_tree_distance()

        # If the tree_distance is a float then
        # the graph is an invalid edge
        # assert self.tree_distance == test_value, 'Tree distance did not match. {} != {}'.format(repr(self.tree_distance),repr(test_value))


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
