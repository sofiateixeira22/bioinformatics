"""
Bioinformatics Course
Group Assignment I
Group E

Members:
- Ana Sofia de Castro Teixeira - 201906031 - MIERSI
- Guilherme Manuel Carvalho de Melo Duarte - 201905583 - M:CC
- José Miguel Ferreira Carvalho - 202005827 - MIERSI
"""

from sys import argv


class DNA:
    def __init__(self, file_path):
        self.file_path = file_path

    def read(self):
        """Reads a FASTA file and returns the sequence."""
        dna_sequence = ''
        with open(self.file_path, 'r') as file:
            for line in file:
                if line[0] != '>':
                    dna_sequence += line.strip()
        file.close()
        return dna_sequence

    def reverse_complement(self, dna_sequence):
        """Returns the reverse complement of a DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(dna_sequence))


class Statistics:
    """Contains the implementations for exercises 1-6."""
    def get(self, dna_sequence):
        """Calculates several statistics from a given DNA sequence."""
        codon_count = {}
        for i in range(len(dna_sequence)-2):
            codon = dna_sequence[i:i+3]
            codon_count[codon] = codon_count.get(codon, 0) + 1

        length = len(dna_sequence)
        frequency = {base: round(dna_sequence.count(base) / length * 100, 2) for base in 'ACGT'}
        gc_content = round((frequency['G'] + frequency['C']) / length * 100, 2)
        start_codons = dna_sequence.count('ATG')
        stop_codons = dna_sequence.count('TAA') + dna_sequence.count('TAG') + dna_sequence.count('TGA')
        most_frequent_codon = max(codon_count, key=codon_count.get)
        least_frequent_codon = min(codon_count, key=codon_count.get)

        return {'length': length,
                'frequency': frequency,
                'gc_content': gc_content,
                'start_codons': start_codons,
                'stop_codons': stop_codons,
                'most_frequent_codon': most_frequent_codon,
                'least_frequent_codon': least_frequent_codon}

    def print(self, stats, strand):
        """Prints given statistics of a DNA sequence."""
        if strand == '+':
            print('Positive Strand Statistics:')
        else:
            print('\nNegative Strand Statistics:')

        print(f'1. Length of the sequence: {stats["length"]}')
        print(f'2. Frequency (in %) of A, C, G, T: {stats["frequency"]}')
        print(f'3. GC content: {stats["gc_content"]}%')
        print(f'4. Number of Start (AUG) codons found: {stats["start_codons"]}')
        print(f'5. Number of Stop Codons (UAA, UAG, UGA): {stats["stop_codons"]}')
        print(f'6. Most and least frequent codons: {stats["most_frequent_codon"]}, {stats["least_frequent_codon"]}')


class Orfs:
    """Contains the implementations for exercises 7-8."""
    def __init__(self):
        self.proteins_file_path = 'all_potential_proteins.txt'
        self.coordinates_file_path = 'orf_coordinates.txt'

    def find(self, dna_sequence, strand):
        """Finds every Open Reading Frame from a DNA sequence, returning a list of protein sequences and correspondent offsets"""
        start_codon = 'ATG'
        stop_codons = ['TAA', 'TAG', 'TGA']
        orfs = []
        protein_sequences = []
        for i in range(len(dna_sequence)-2):
            if dna_sequence[i:i+3] != start_codon:
                continue
            for j in range(i, len(dna_sequence)-2, 3):
                if dna_sequence[j:j+3] not in stop_codons:
                    continue
                orf_length = j + 3 - i
                if orf_length >= 150:
                    orfs.append((i+1, j+3))  # Adding 1 to start for 1-based indexing
                    protein_sequences.append(self._translation(dna_sequence[i:j+3]))
                break
        if strand == '-': # Is this a thing?
            orfs = [(len(dna_sequence) - orf[1] + 1, len(dna_sequence) - orf[0] + 1) for orf in orfs]
        return orfs, protein_sequences

    def _translation(self, dna_sequence):
        """Translates sequence to protein using an internal dictionary with the standard genetic code."""
        tc = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'TGT': 'C', 'TGC': 'C',
              'GAT': 'D', 'GAC': 'D',
              'GAA': 'E', 'GAG': 'E',
              'TTT': 'F', 'TTC': 'F',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
              'CAT': 'H', 'CAC': 'H',
              'ATA': 'I', 'ATT': 'I', 'ATC': 'I',
              'AAA': 'K', 'AAG': 'K',
              'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'ATG': 'M', 'AAT': 'N', 'AAC': 'N',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'TGG': 'W',
              'TAT': 'Y', 'TAC': 'Y',
              'TAA': '_', 'TAG': '_', 'TGA': '_'}
        protein_sequence = ''
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in tc:
                protein_sequence += tc[codon]
            else:
                raise Exception("Invalid DNA sequence.")  # Right?
        return protein_sequence

    def read(self):
        """Reads from file the coordinates of the protein sequences found."""
        orfs = []
        with open(self.coordinates_file_path, 'r') as orf_file:
            for line in orf_file:
                line = line.strip().split(',')
                orf_start = int(line[0])
                orf_end = int(line[1])
                orfs.append((orf_start, orf_end))
        return orfs

    def save(self, orfs, protein_sequences, strand):
        """Stores Protein sequences and respective coordinates from a DNA sequence into files."""
        with open(self.proteins_file_path, 'w+') as proteins_file, open(self.coordinates_file_path, 'w+') as coordinates_file:
            for i, (orf, protein) in enumerate(zip(orfs, protein_sequences), start=1):
                proteins_file.write(f'{protein}\n')
                coordinates_file.write(f'{orf[0]}, {orf[1]}, ORF{i}\n')
        proteins_file.close()
        coordinates_file.close()
        print(f'\nIdentified and saved protein sequences for {"positive" if strand == "+" else "negative"} ORFs.')


class Overlap:
    """Contains the implementations for exercise 9."""
    def __init__(self, file_path):
        self.file_path = file_path
    def read_exons(self):
        """Extracts every id, start and end position of an exon from an annotation file."""
        exons = []
        with open(self.file_path, 'r') as file:
            for line in file:
                line = line.split('\t')
                type = line[2]
                if type != 'exon':
                    continue
                start_codon = int(line[3])
                end_codon = int(line[4])
                gene_id = line[8].split(';')[0].split(' ')[1].replace('"', '')
                exons.append((gene_id, start_codon, end_codon))
        return exons

    def _calculate_overlap(self, orfA_start, orfA_end, orfB_start, orfB_end):
        """Computes the overlap between two sets of indexes."""
        start = max(orfA_start, orfB_start)
        end = min(orfA_end, orfB_end)
        intersection = max(end - start + 1, 0)
        overlap = round(intersection / (orfA_end - orfA_start + 1) * 100, 2)
        return overlap

    def get(self):
        """For every ORF in an annotation file, computes the longest overlap with the ORFs found."""
        exons_annotation = self.read_exons()
        exons_found = Orfs().read()
        overlap = []
        for (orfA_id, orfA_start, orfA_end) in exons_annotation:
            longest_overlap = 0
            for (orfB_start, orfB_end) in exons_found:
                overlap_AB = self._calculate_overlap(orfA_start, orfA_end, orfB_start, orfB_end)
                longest_overlap = max(longest_overlap, overlap_AB)
            overlap.append((orfA_id, longest_overlap))
        return overlap

    def print(self, overlap):
        """Prints the longest overlap between the ORFs in an annotation file, and the ones found on a DNA sequence."""
        print(f'\nLongest overlap for each ORF in the annotation file {self.file_path}:')
        for (orf_id, longest_overlap) in overlap:
            print(f'\t{orf_id}:\t{longest_overlap}%')


if __name__ == '__main__':
    # Get both strands from a DNA sequence
    fasta_file = argv[1]
    dna = DNA(fasta_file)
    dna_sequence = dna.read()
    reverse_dna_sequence = dna.reverse_complement(dna_sequence)

    # Get statistics for both strands
    statistics = Statistics()
    positive_statistics = statistics.get(dna_sequence)
    negative_statistics = statistics.get(reverse_dna_sequence)
    statistics.print(positive_statistics, strand='+')
    statistics.print(negative_statistics, strand='-')

    # Identify ORFs for both strands
    orfs = Orfs()
    positive_orfs, positive_protein_sequences = orfs.find(dna_sequence, strand='+')
    negative_orfs, negative_protein_sequences = orfs.find(reverse_dna_sequence, strand='-')
    orfs.save(positive_orfs, positive_protein_sequences, strand='+')
    orfs.save(negative_orfs, negative_protein_sequences, strand='-')

    # Overlap with annotation
    annotation_file = argv[2]
    overlap = Overlap(annotation_file)
    annotation_overlap = overlap.get()
    overlap.print(annotation_overlap)
