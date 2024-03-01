def read_fasta(file_path):
    """Reads a FASTA file and returns the sequence."""
    sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    file.close()
    return sequence

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def calculate_statistics(seq):
    """Calculates various statistics from a given DNA sequence."""
    length = len(seq)
    freq = {base: round(seq.count(base) / length * 100, 2) for base in 'ACGT'}
    gc_content = round((freq['G'] + freq['C']) / 100, 2)
    start_codons = seq.count('ATG')
    stop_codons = seq.count('TAA') + seq.count('TAG') + seq.count('TGA')
    
    codon_count = {}
    for i in range(0, len(seq), 3):
        # checks if codon has 3 bases (end of sequence codons may not have 3 bases)
        if(len(seq[i:i+3]) == 3):
           codon = seq[i:i+3]
           codon_count[codon] = codon_count.get(codon, 0) + 1

    most_frequent_codon = max(codon_count, key=codon_count.get), codon_count[max(codon_count, key=codon_count.get)]
    least_frequent_codon = min(codon_count, key=codon_count.get), codon_count[min(codon_count, key=codon_count.get)]
    
    return length, freq, gc_content, start_codons, stop_codons, most_frequent_codon, least_frequent_codon

def find_orfs(sequence):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    protein_sequences = []
    for i in range(len(sequence)):
        if sequence[i:i+3] == start_codon:
            for j in range(i, len(sequence), 3):
                if sequence[j:j+3] in stop_codons:
                    orf_length = j + 3 - i
                    if orf_length >= 150:
                        orfs.append((i+1, j+3))  # Adding 1 to start for 1-based indexing
                        protein_sequences.append(translate_to_protein(sequence[i:j+3]))
                    break
    return orfs, protein_sequences

def translate_to_protein(dna_sequence):
    # Translates sequence to protein using an internal dictionary with the standard genetic code.
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    
    translated = ''
    for i in range(0, len(dna_sequence), 3):
        if(dna_sequence[i:i+3] in tc):
            translated += tc[dna_sequence[i:i+3]]

    return translated

def save_orf_information(orfs, protein_sequences):
    with open('all_potential_proteins.fasta', 'w') as proteins_file, open('orf_coordinates.txt', 'w') as coordinates_file:
        for idx, (orf, protein) in enumerate(zip(orfs, protein_sequences), start=1):
            proteins_file.write(f'>Protein_ORF{idx}\n{protein}\n')
            coordinates_file.write(f'{orf[0]}, {orf[1]}, ORF{idx}\n')
    proteins_file.close()
    coordinates_file.close()

def calculate_overlap(start, end):
    max_overlap = 0
    with open('orf_coordinates.txt', 'r') as orf_file:
        for line in orf_file:
            line = line.strip().split(',')
            orf_start = int(line[0])
            orf_end = int(line[1])
            if(orf_start >= start and orf_end <= end):
                overlap = ((orf_end-orf_start) / (end-start)) * 100
                max_overlap = max(max_overlap, overlap)
    return max_overlap

# Main script
fasta_file = 'sequence_chr1.fasta'
sequence = read_fasta(fasta_file)
reverse_sequence = reverse_complement(sequence)

# Calculate statistics for both strands
stats_positive = calculate_statistics(sequence)
stats_negative = calculate_statistics(reverse_sequence)

# Output statistics
print("Positive Strand Statistics:")
print(f"1. Length of the sequence: {stats_positive[0]}")
print(f"2. Frequency (in %) of A, C, G, T: {stats_positive[1]}")
print(f"3. GC content: {stats_positive[2]}")
print(f"4. Number of Start (AUG) codons found: {stats_positive[3]}")
print(f"5. Number of Stop Codons (UAA, UAG, UGA): {stats_positive[4]}")
print(f"6. Most and least frequent codons: {stats_positive[5][0]} {stats_positive[5][1]}, {stats_positive[6][0]}, {stats_positive[6][1]}")

print("\nNegative Strand Statistics:")
print(f"1. Length of the sequence: {stats_negative[0]}")
print(f"2. Frequency (in %) of A, C, G, T: {stats_negative[1]}")
print(f"3. GC content: {stats_negative[2]}")
print(f"4. Number of Start (AUG) codons found: {stats_negative[3]}")
print(f"5. Number of Stop Codons (UAA, UAG, UGA): {stats_negative[4]}")
print(f"6. Most and least frequent codons: {stats_negative[5][0]} {stats_negative[5][1]}, {stats_negative[6][0]} {stats_negative[6][1]}")

# Identify ORFs for both strands
orfs_positive, protein_sequences_positive = find_orfs(sequence)
orfs_negative, protein_sequences_negative = find_orfs(reverse_sequence)

# Adjust coordinates for negative strand ORFs
orfs_negative_adjusted = [(len(sequence) - orf[1] + 1, len(sequence) - orf[0] + 1) for orf in orfs_negative]

# Save ORF information
save_orf_information(orfs_positive, protein_sequences_positive)
save_orf_information(orfs_negative_adjusted, protein_sequences_negative)

print(f"\nIdentified and saved information for {len(orfs_positive) + len(orfs_negative)} ORFs.")

# Overlap with annotation
with open('genes_chr1.gtf', 'r') as genes_file:
    for line in genes_file:
        line = line.strip().split('\t')
        type = line[2]
        if(type == 'exon'):
            start_codon = int(line[3])
            end_codon = int(line[4])
            percentage = calculate_overlap(start_codon, end_codon)
            gene_id = line[8].split(';')[0].split(' ')[1].replace('"', '')
            print(f"{gene_id} {percentage}%")