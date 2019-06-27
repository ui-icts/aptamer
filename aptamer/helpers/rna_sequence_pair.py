import RNA
import Levenshtein


class RNASequencePair(object):
    """Graph edge. Pair of RNASequence objects."""
    def __init__(self, seq1, seq2):
        self.sequence1 = seq1
        self.sequence2 = seq2
        self.energy_delta = None
        self.edit_distance = None
        self.tree_distance = None
        self.is_valid_edge = False  # used only with seed option

        self.update_energy_delta()
        self.update_edit_distance()
        self.update_tree_distance()

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

        test_tree_distance = RNA.tree_edit_distance(
                RNA.make_tree(RNA.expand_Full(self.sequence1.structure)),
                RNA.make_tree(RNA.expand_Full(self.sequence2.structure))
                                                )

        self.tree_distance = test_tree_distance

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
