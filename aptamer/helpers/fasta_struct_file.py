import re
import itertools
from helpers.rna_sequence import RNASequence


class FastaStructFile(object):

    def __init__(self, in_fh):
        self.in_fh = in_fh
        self.cluster_size_re = re.compile('SIZE=(\d+)')

    def rna_seq_objs(self, calc_structures, prefix, suffix):
        while True:
            try:
                # need to move through a triplet file structure, not fasta
                header, sequence, structure = list(
                    itertools.islice(self.in_fh, 3)
                )
            except ValueError:  # end of file
                break
            sequence = sequence.strip('\n\r')
            structure = structure.strip('\n\r')

            if not structure.count('(') == structure.count(')'):
                continue

            if calc_structures:
                sequence = '%s%s%s'.replace('T', 'U') % (
                    prefix, sequence, suffix
                )
            else:
                sequence = sequence.replace('T', 'U')

            try:
                cluster_size = self.cluster_size_re.search(header)
                cluster_size = cluster_size.group(1)
            except AttributeError:
                print('Not able to find cluster size. Setting to 1.')

            if cluster_size is None:
                cluster_size = 1

            header = header.replace('>', '').strip()
            curr_seq = RNASequence(header, cluster_size, sequence)
            curr_seq.free_energy = 1
            curr_seq.ensemble_free_energy = 1
            curr_seq.ensemble_probability = 1
            curr_seq.ensemble_diversity = 1
            curr_seq.structure = structure

            yield curr_seq

    def write_combinations(self, output_file_name):
        objs = self.rna_seq_objs(False,None,None)
        seq_pairs = (p for p in itertools.combinations(objs, 2))

        with open(output_file_name, 'w') as out:
            for p1, p2 in seq_pairs:
                out.write(p1.structure + '\n')
                out.write(p2.structure + '\n')

