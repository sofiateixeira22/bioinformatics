def read_fasta(file_path):
    """Reads a FASTA file and returns the sequence."""
    sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
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
        codon = seq[i:i+3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            codon_count[codon] = 1

    most_frequent_codon = max(codon_count, key=codon_count.get), max(codon_count.values())
    least_frequent_codon = min(codon_count, key=codon_count.get), min(codon_count.values())
    
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
    # Simple translation (mock function for demonstration; doesn't cover real translation)
    return ''.join(['X' for _ in range(0, len(dna_sequence), 3)])

def save_orf_information(orfs, protein_sequences):
    with open('all_potential_proteins.fasta', 'w') as proteins_file, open('orf_coordinates.txt', 'w') as coordinates_file:
        for idx, (orf, protein) in enumerate(zip(orfs, protein_sequences), start=1):
            proteins_file.write(f'>Protein_ORF{idx}\n{protein}\n')
            coordinates_file.write(f'{orf[0]}, {orf[1]}, ORF{idx}\n')

# Main script
fasta_file = 'BioInformatica/bioinformatics/sequence_chr1.fasta'
sequence = read_fasta(fasta_file)
reverse_sequence = reverse_complement(sequence)

# Calculate statistics for both strands
stats_positive = calculate_statistics(sequence)
stats_negative = calculate_statistics(reverse_sequence.replace('T', 'U'))

# Output statistics
print("Positive Strand Statistics:")
print(f"1. Length of the sequence: {stats_positive[0]}")
print(f"2. Frequency (in %) of A, C, G, T: {stats_positive[1]}")
print(f"3. GC content: {stats_positive[2]}")
print(f"4. Number of Start (AUG) codons found: {stats_positive[3]}")
print(f"5. Number of Stop Codons (UAA, UAG, UGA): {stats_positive[4]}")
print(f"6. Most and least frequent codons: {stats_positive[5][0]} ({stats_positive[5][1]}), {stats_positive[6][0]} ({stats_positive[6][1]})")

print("\nNegative Strand Statistics:")
print(f"1. Length of the sequence: {stats_negative[0]}")
print(f"2. Frequency (in %) of A, C, G, U: {stats_negative[1]}")  # Note the U for RNA
print(f"3. GC content: {stats_negative[2]}")
print(f"4. Number of Start (AUG) codons found: {stats_negative[3]}")
print(f"5. Number of Stop Codons (UAA, UAG, UGA): {stats_negative[4]}")
print(f"6. Most and least frequent codons: {stats_negative[5][0]} ({stats_negative[5][1]}), {stats_negative[6][0]} ({stats_negative[6][1]})")

# Identify ORFs for both strands
orfs_positive, protein_sequences_positive = find_orfs(sequence)
orfs_negative, protein_sequences_negative = find_orfs(reverse_sequence)

# Adjust coordinates for negative strand ORFs
orfs_negative_adjusted = [(len(sequence) - orf[1] + 1, len(sequence) - orf[0] + 1) for orf in orfs_negative]

# Save ORF information
save_orf_information(orfs_positive, protein_sequences_positive)
save_orf_information(orfs_negative_adjusted, protein_sequences_negative)

print(f"\nIdentified and saved information for {len(orfs_positive) + len(orfs_negative)} ORFs.")
