import re
import itertools
import subprocess
import concurrent.futures
from Bio import SeqIO
from helpers.rna_sequence import RNASequence

class FastaFile(object):

    def __init__(self, in_fh):
        self.in_fh = in_fh
        self.cluster_size_re = re.compile('SIZE=(\d+)')
        self.prefix = ''
        self.suffix = ''
        self.run_mfold = False
        self.pass_options = ''
        self.vienna_version = 2

    def rna_seq_objs(self):
        """Process input file as fasta. Populate RNASequence (graph vertex)
        objects.
        """

        in_fh = self.in_fh
        cluster_size_re = self.cluster_size_re
        prefix = self.prefix
        suffix = self.suffix
        run_mfold = self.run_mfold
        pass_options = self.pass_options
        vienna_version = self.vienna_version

        folds = []
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for record in SeqIO.parse(in_fh, 'fasta'):
                sequence = '%s%s%s'.replace('T', 'U') % (
                    prefix, str(record.seq), suffix
                )
                cluster_size = 1
                try:
                    cluster_size = cluster_size_re.search(record.description)
                    cluster_size = cluster_size.group(1)
                except AttributeError:
                    pass
                    # print('Not able to find cluster size. Setting to 1.')

                if cluster_size is None:
                    cluster_size = 1

                # find structure
                curr_seq = RNASequence(record.id, cluster_size, sequence)
                if run_mfold:
                    mfold_f = executor.submit(run_mfold_prg, sequence, pass_options)
                    folds.append(mfold_f)
                else:
                    rnafold_f = executor.submit(run_rnafold, sequence, vienna_version, pass_options)
                    folds.append(rnafold_f)

        for folded_output in folds:
            if run_mfold:
                curr_seq.structure, curr_seq.energy_dict = folded_output.result()
                curr_seq.free_energy = curr_seq.energy_dict['dG']
            else:
                rnafold_out = folded_output.result()
                rnafold_out = rnafold_out.split('\n')
                try:
                    curr_seq.structure, curr_seq.free_energy = (
                        rnafold_out[1].split(' (')
                    )
                except (ValueError, IndexError):
                    print('Error running RNAfold:\n%s\nExiting.' % rnafold_out)
                    sys.exit(1)

                print('%s\n' % rnafold_out)
                try:
                    curr_seq.free_energy = abs(
                        float(curr_seq.free_energy.replace(')', ''))
                    )
                    curr_seq.ensemble_free_energy = abs(
                        float(rnafold_out[2].split('[')[1].replace(']', ''))
                    )
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
                except IndexError:
                    print(
                        'Error parsing RNAfold output. '
                        '(Couldn\'t find statistics.) Please check '
                        'RNAfold options.'
                    )
                    sys.exit(1)

            (yield curr_seq)


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
    return stdout_value


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
    return structure, energy_stats
