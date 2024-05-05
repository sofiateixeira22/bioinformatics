import sys
import requests
import pprint  # to remove later - pretty print for json objects

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/search?query='


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

    def _save_sequence_fasta(self):
        print('> Saving sequence in current directory as sequence.fasta')
        uniProt_id = self.seq_info['uniProtkbId']
        sequence = self.seq_info['sequence']['value']
        with open('sequence.fasta', 'w') as f:
            f.write('>' + uniProt_id + '\n')
            f.write(sequence)

    def fetch_sequence_information(self):
        seq_info = self._fetch_UniProt()
        if not seq_info:
            print(f'> No results found with the sequence ID {seq_id}')
        else:
            self.seq_info = seq_info
            uniProt_id = self.seq_info['uniProtkbId']
            print(f'> Found match on UniProt database with ID {uniProt_id}')
            self._save_sequence_fasta()


if __name__ == '__main__':
    # seq_id = sys.argv[1]
    seq_id = 'P68871'
    phylogenetics = Phylogenetics(seq_id)
    phylogenetics.fetch_sequence_information()
