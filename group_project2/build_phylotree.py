import sys
import requests
import pprint  # to remove later - pretty print for json objects

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
# from Bio import Entrez

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/search?query='
MAX_SPECIES = 10
SEQUENCE_FASTA_FILE = 'sequence.fasta'
TOP_SEQUENCE_FASTA_FILE = 'sequences_to_analyse.fasta'

# def print_sequence_info(seq_info):
#     print(f'--- Sequence Info ---')
#     uniprot_id = seq_info['results'][0]['uniProtkbId']
#     print(f'UniProtKB ID: {uniprot_id}')
#     first_pubic_date = seq_info['results'][0]['entryAudit']['firstPublicDate']
#     print(f'First Public Date: {first_pubic_date}')
#     organism = seq_info['results'][0]['organism']['scientificName']
#     print(f'Organism: {organism}')
#     recommended_name = seq_info['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
#     print(f'Recommended Name: {recommended_name}')
#     gene_name = seq_info['results'][0]['genes'][0]['geneName']['value']
#    print(f'Gene Name: {gene_name}')

class Phylogenetics:
    def __init__(self, seq_id):
        self.seq_id = seq_id
        self.seq_info = None

    def _fetch_UniProt(self):
        try:
            query = UNIPROT_URL + self.seq_id
            response = requests.get(query)
            results = response.json()
            sequence = results['results'][0]
            return sequence
        except:
            return None

    @staticmethod
    def _save_sequence_fasta(seq_ids, sequences, output_file):
        print(f'\t> Saving sequence{"s" if len(seq_ids) > 1 else ""} in {output_file}')
        with open(output_file, 'w') as f:
            for seq_id, sequence in zip(seq_ids, sequences):
                f.write('>' + seq_id + '\n')
                f.write(sequence + '\n')

    def data_collection(self):
        seq_info = self._fetch_UniProt()
        if not seq_info:
            print(f'> No results found with the sequence ID {seq_id}')
        else:
            self.seq_info = seq_info
            uniProt_id = self.seq_info['uniProtkbId']
            print(f'> Found match on UniProt database with ID {uniProt_id}')
            sequence = self.seq_info['sequence']['value']
            self._save_sequence_fasta([uniProt_id], [sequence], SEQUENCE_FASTA_FILE)

    def BLAST_analysis(self):
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
        self._save_sequence_fasta(best_species, best_species_sequence, TOP_SEQUENCE_FASTA_FILE)


if __name__ == '__main__':
    # seq_id = sys.argv[1]
    seq_id = 'P68871'

    phylogenetics = Phylogenetics(seq_id)
    phylogenetics.data_collection()
    phylogenetics.BLAST_analysis()
