import sys
import requests
from json import loads

def fetch_seq_info(seq_id):
    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + seq_id
    response = requests.get(url)
    seq_info = response.json()
    sequence = seq_info['results'][0]['sequence']['value']
    print_sequence_info(seq_info)
    write_sequence_fasta(seq_id, sequence)

def print_sequence_info(seq_info):
    print(f'--- Sequence Info ---')
    uniprot_id = seq_info['results'][0]['uniProtkbId']
    print(f'UniProtKB ID: {uniprot_id}')
    first_pubic_date = seq_info['results'][0]['entryAudit']['firstPublicDate']
    print(f'First Public Date: {first_pubic_date}')
    organism = seq_info['results'][0]['organism']['scientificName']
    print(f'Organism: {organism}')
    recommended_name = seq_info['results'][0]['proteinDescription']['recommendedName']['fullName']['value']
    print(f'Recommended Name: {recommended_name}')
    gene_name = seq_info['results'][0]['genes'][0]['geneName']['value']
    print(f'Gene Name: {gene_name}')
    
def write_sequence_fasta(seq_identifier, sequence):
    with open('sequence.fasta', 'w') as f:
        f.write('>' + seq_identifier + '\n')
        f.write(sequence)

if __name__ == '__main__':
    # seq_identifier = sys.argv[1]
    seq_identifier = 'P68871'
    fetch_seq_info(seq_identifier)