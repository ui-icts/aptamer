class RNASequence:
    __slots__ = 'name', 'cluster_size', 'sequence', 'structure', 'free_energy', 'energy_dict', 'ensemble_free_energy', 'ensemble_probability', 'ensemble_diversity', 'use_for_comparison'

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

    def output(self):
        print(">%s  " % (self.name))
        print(self.sequence)
        print(self.structure)

    def __str__(self):
        return '\n'.join(
            ['%s : %s' % (z, getattr(self, z)) for z in self.__slots__]
        )

    __repr__ = __str__
