import re
import itertools
import subprocess
import concurrent.futures
from multiprocessing import Pool
from Bio import SeqIO
from helpers.rna_sequence import RNASequence

class FastaFile(object):

    def __init__(self, input_file_name):
        self.input_file_name = input_file_name
        self.cluster_size_re = re.compile('SIZE=(\d+)')
        self.prefix = ''
        self.suffix = ''
        self.run_mfold = False
        self.pass_options = ''
        self.vienna_version = 2

    def __iter__(self):
        """Process input file as fasta. Populate RNASequence (graph vertex)
        objects.
        """

        with open(self.input_file_name, 'r') as in_fh:
            for record in SeqIO.parse(in_fh, 'fasta'):
                seq = self.build_rna_sequence(record)
                (yield seq)


    def p_seq_objs(self):
        pool = Pool()
        with open(self.input_file_name, 'r') as in_fh:
            records = SeqIO.parse(in_fh, 'fasta')
            return pool.map(self.build_rna_sequence, records)


    def get_cluster_size(self, record):
        cluster_size_re = self.cluster_size_re
        cluster_size = 1
        try:
            cluster_size = cluster_size_re.search(record.description)
            cluster_size = cluster_size.group(1)
        except AttributeError:
            pass
            # print('Not able to find cluster size. Setting to 1.')

        if cluster_size is None:
            cluster_size = 1

        return cluster_size

    def get_sequence(self, record):
        prefix = self.prefix
        suffix = self.suffix
        sequence = '%s%s%s'.replace('T', 'U') % (
            prefix, str(record.seq), suffix
        )
        return sequence

    def build_rna_sequence(self, record):
        sequence = self.get_sequence(record)
        cluster_size = self.get_cluster_size(record)

        # find structure
        curr_seq = RNASequence(record.id, cluster_size, sequence)

        fold_output = self.mfold(sequence) if self.run_mfold else self.rnafold(sequence)
        fold_output.update(curr_seq)

        return curr_seq

    def rnafold(self, sequence):
        return run_rnafold(sequence, self.vienna_version, self.pass_options)

    def mfold(self, sequence):
        return run_mfold_prg(sequence,self.pass_options)

class RNAFoldOutput(object):

    def __init__(self, rnafold_out):
        rnafold_out = rnafold_out.split('\n')
        try:
            structure, free_energy = (
                rnafold_out[1].split(' (')
            )
            self.structure = structure
            self.free_energy = abs( float(free_energy.replace(')', '')))
        except (ValueError, IndexError):
            print('Error running RNAfold:\n%s\nExiting.' % rnafold_out)
            sys.exit(1)

        print('%s\n' % rnafold_out)

        try:
            self.ensemble_free_energy = abs(
                float(rnafold_out[2].split('[')[1].replace(']', ''))
            )
            self.ensemble_probability = abs(float(
                rnafold_out[4].split(';')[0].replace(
                    ' frequency of mfe structure in ensemble ', ''
                )
            ))
            self.ensemble_diversity = abs(float(
                rnafold_out[4].split(';')[1].replace(
                    ' ensemble diversity ', ''
                )
            ))
        except IndexError:
            print(
                'Error parsing RNAfold output. '
                '(Couldn\'t find statistics.) Please check '
                'RNAfold options.'
            )
            sys.exit(1)

    def update(self, rna_sequence):
        rna_sequence.structure = self.structure
        rna_sequence.free_energy = self.free_energy
        rna_sequence.free_energy = self.free_energy
        rna_sequence.ensemble_free_energy = self.ensemble_free_energy
        rna_sequence.ensemble_probability = self.ensemble_probability
        rna_sequence.ensemble_diversity = self.ensemble_diversity

def run_rnafold(seq, vienna_version, pass_options):
    print('##################')
    print('Running RNAFold...')
    print('##################')
    cmd = None
    if pass_options is not None:
        cmd = ['RNAfold %s' % pass_options]
    elif vienna_version == 1:
        cmd = ['RNAfold -p -T 30 -noLP -noPS -noGU']
    elif vienna_version == 2:
        cmd = ['RNAfold -p -T 30 --noLP --noPS --noGU']
    sffproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        stdin=subprocess.PIPE, close_fds=True, shell=True,
        encoding='ascii'
    )
    print('Command:', cmd[0])
    stdout_value, stderr_value = sffproc.communicate(seq)
    return RNAFoldOutput(stdout_value)


class MFoldOutput(object):
    def __init__(self, mfold_out):
        self.structure, self.energy_dict = mfold_out
        self.free_energy = curr_seq.energy_dict['dG']

    def update(self, rna_sequence):
        rna_sequence.structure = self.structure
        rna_sequence.energy_dict = self.energy_dict
        rna_sequence.free_energy = self.energy_dict['dG']


def run_mfold_prg(seq, pass_options):
    print('##################')
    print('Running mfold...')
    print('##################')
    if not os.path.exists('mfold_out'):
        os.mkdir('mfold_out')
    os.chdir('mfold_out')
    temp_filename = 'mfold_temp.txt'
    with open(temp_filename, 'w') as f:
        f.write('%s\n' % seq)

    if pass_options is not None:
        cmd_string = 'mfold SEQ=%s %s' % (temp_filename, pass_options)
    else:
        cmd_string = 'mfold SEQ=%s T=30' % temp_filename

    ret = subprocess.call(
        cmd_string, stderr=subprocess.STDOUT, shell=True
    )

    if ret != 0:
        print(
            'Error when running mfold. Return code %d. '
            'See mfold log file for details.' % ret
        )
        sys.exit(ret)
    print()
    structure = convert_ct_to_bracket_dot('%s.ct' % temp_filename)
    energy_stats = get_mfold_stats('%s.det' % temp_filename)
    os.chdir('..')

    return MFoldOutput((structure, energy_stats))


def convert_ct_to_bracket_dot(ct_filename):
    bracket_dot = ''
    with open(ct_filename, 'r') as f:
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


def get_mfold_stats(det_filename):
    # fixed values within line
    valid_pattern = {
        0: 'dG', 1: '=', 3: 'dH', 4: '=', 6: 'dS', 7: '=', 9: 'Tm', 10: '='
    }
    energy_stats = {}
    for key, value in valid_pattern.items():
        if value != '=':
            energy_stats[value] = None
    with open(det_filename, 'r') as f:
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


