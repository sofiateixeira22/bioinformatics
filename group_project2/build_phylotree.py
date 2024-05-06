from requests import get, post
from sys import argv
from time import sleep

import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO
from Bio import Phylo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/'
EBI_URL = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'

SEQUENCE_FASTA_FILE = 'sequence.fasta'
TOP_SEQUENCE_FASTA_FILE = 'sequences_to_analyse.fasta'
ALIGNMENT_FILE = 'alignment.txt'
TREE_FILE = 'phylotree.png'

MAX_SPECIES = 10


def do_request(method, url, command=None, parameters=None):
    response = None
    if method == 'get':
        response = get(url+command, data=parameters)
    elif method == 'post':
        response = post(url + command, data=parameters)

    if response.status_code == 200:
        return response
    else:
        print(f'\t> Request failed with status code {response.status_code}')
        print(f'\t> Response:\n{response.text}\n', response.text)
        return None


def save_fasta(seq_ids, sequences, output_file):
    print(f'\t> Saving sequence{"s" if len(seq_ids) > 1 else ""} in {output_file}')
    with open(output_file, 'w') as f:
        for seq_id, sequence in zip(seq_ids, sequences):
            f.write(f'>{seq_id.replace(" ", "-")}\n{sequence}\n')


class Phylogenetics:
    def __init__(self, seq_id):
        self.seq_id = seq_id

    def data_collection(self):
        try:
            results = do_request('get', UNIPROT_URL, f'/search?query={self.seq_id}')
            seq_info = results.json()['results'][0]
            species = seq_info['organism']['scientificName']
            print(f'> Found match on UniProt database for species {species}')
            sequence = seq_info['sequence']['value']
            save_fasta([species], [sequence], SEQUENCE_FASTA_FILE)
        except:
            print(f'> No results found with the sequence ID {seq_id}')

    @staticmethod
    def BLAST_analysis():
        print('> Performing a BLAST search')
        sequence = SeqIO.read(SEQUENCE_FASTA_FILE, 'fasta').seq

        # To delete, but so we understand:
        #   - blastp: protein blast, for protein vs protein comparison
        #   - nr: default protein database, stands for non-redundant protein sequences

        # print('\t> Queueing server')
        # result_handle = NCBIWWW.qblast(program='blastp', database='nr', sequence=sequence, hitlist_size=250)
        # with open('blast_results.xml', 'w') as out_handle:
        #     out_handle.write(result_handle.read())

        print('\t> Parsing result')
        with open('blast_results.xml', 'r') as result_handle:
            blast_record = NCBIXML.read(result_handle)

        species_score = {}
        species_sequence = {}
        for alignment in blast_record.alignments:
            species = alignment.hit_def.split('[')[1].split(']')[0]
            if species in ['Homo sapiens', 'synthetic construct']:
                continue

            for hsp in alignment.hsps:
                sbj_score = hsp.score
                sbj_sequence = hsp.sbjct

                if (species in species_score and species_score[species] < sbj_score) or species not in species_score:
                    species_score[species] = sbj_score
                    species_sequence[species] = sbj_sequence

        best_species = sorted(species_score.items(), key=lambda x: x[1], reverse=True)
        best_species = [species[0] for species in best_species[:(min(len(best_species), MAX_SPECIES))]]
        best_species_sequence = [species_sequence[species] for species in best_species]
        save_fasta(best_species, best_species_sequence, TOP_SEQUENCE_FASTA_FILE)

    @staticmethod
    def multiple_sequence_alignment():
        print('> Applying Multiple Sequence Alignment')
        sequence = SeqIO.read(SEQUENCE_FASTA_FILE, 'fasta')
        species = list(SeqIO.parse(TOP_SEQUENCE_FASTA_FILE, 'fasta'))

        run_response = do_request('post', EBI_URL, '/run', parameters={
            'email': 'abc@abc.com',
            'sequence': '\n'.join(f'>{seq_rec.name}\n{seq_rec.seq}' for seq_rec in [sequence] + species)
        })
        if not run_response:
            return
        job_id = run_response.text
        print(f'\t> Job successfully submitted with ID {job_id}')

        print(f'\t> Waiting for result...')
        while True:
            status_response = do_request('get', EBI_URL, f'/status/{job_id}')
            if not status_response:
                return
            status = status_response.text
            if status == 'FINISHED':
                break
            sleep(1)
        print(f'\t> Job done')

        result_response = do_request('get', EBI_URL, f'/result/{job_id}/aln-clustal_num')
        if not result_response:
            return
        result = result_response.text

        with open(ALIGNMENT_FILE, "w") as f:
            f.write(result)
        print(f'\t> Saving results on {ALIGNMENT_FILE}')

    @staticmethod
    def phylogenetic_tree():
        print('> Building Phylogenetic Tree')
        alignment = AlignIO.read(ALIGNMENT_FILE, 'clustal')

        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        tree_constructor = DistanceTreeConstructor()
        tree = tree_constructor.upgma(distance_matrix)

        plt.rcParams.update({'font.size': 11})
        figure = plt.figure(figsize=(15, 8))
        axes = figure.add_subplot(111, xticks=[], yticks=[])
        Phylo.draw(tree, axes=axes, branch_labels=lambda c: round(c.branch_length, 5), do_show=False)
        plt.xlabel('')
        plt.ylabel('')
        plt.savefig(TREE_FILE, dpi=200, bbox_inches='tight')
        print(f'\t> Saving tree on {TREE_FILE}')


if __name__ == '__main__':
    phylogenetics = Phylogenetics(argv[1])
    phylogenetics.data_collection()
    phylogenetics.BLAST_analysis()
    phylogenetics.multiple_sequence_alignment()
    phylogenetics.phylogenetic_tree()
