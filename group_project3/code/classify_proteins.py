from sys import argv, exit
import argparse
from itertools import product

from Bio import SeqIO
import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import make_scorer, precision_score, recall_score, f1_score


def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, 'fasta'):
        sequences[record.id.split('|')[1]] = str(record.seq)
    return sequences


def generate_mers(k):
    aminoacids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return [''.join(pair) for pair in product(aminoacids, repeat=k)]


def compute_ffp(sequences, k_mers):
    ffp_data = []
    number_k_mers = len(k_mers)

    for seq_id, sequence in sequences.items():
        seq_length = len(sequence)
        ffp = {k: 0 for k in k_mers}
        for i in range(seq_length - 1):
            nucleotides = sequence[i:i+len(k_mers[0])]
            if nucleotides in ffp:
                ffp[nucleotides] += 1

        for k in ffp:
            ffp[k] /= number_k_mers

        ffp_data.append([seq_id] + [ffp[k] for k in k_mers])

    columns = ['SequenceID'] + k_mers
    return pd.DataFrame(ffp_data, columns=columns)


def classification(ffp):
    x = ffp.drop(columns=['SequenceID', 'Class'])
    y = ffp['Class']

    classifiers = {
        'SVM': SVC(),
        'Random Forest': RandomForestClassifier(),
        'Naive Bayes': GaussianNB()
    }

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    results = {}
    for name, clf in classifiers.items():
        f1_score = cross_val_score(clf, x, y, cv=cv, scoring='f1')
        results[name] = {
            'precision': cross_val_score(clf, x, y, cv=cv, scoring='precision').mean(),
            'recall': cross_val_score(clf, x, y, cv=cv, scoring='recall').mean(),
            'f1_mean': f1_score.mean(),
            'f1_std': f1_score.std()
        }

    return pd.DataFrame(results).T


if __name__ == '__main__':
    # Argument Parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--file_a', type=str, required=True)
    parser.add_argument('-b', '--file_b', type=str, required=True)
    parser.add_argument('-k', '--k_value', type=int, required=True)

    args = parser.parse_args()
    file_a = args.file_a
    file_b = args.file_b
    k = args.k_value

    # Task 1 - read sequences from fasta files and return data as dictionary with key as sequence ID and value as the sequence
    globin_sequences = read_fasta(file_a)
    zincfinger_sequences = read_fasta(file_b)

    # Task 2 - generate al the k-mers of amino-acids
    k_mers = generate_mers(k)

    # Task 3 - create and fill a pandas dataframe
    globin_ffp = compute_ffp(globin_sequences, k_mers)
    zincfinger_ffp = compute_ffp(zincfinger_sequences, k_mers)

    # Task 4 - add extra column (class) to the table to include the label of the type of protein
    globin_ffp['Class'] = 0
    zincfinger_ffp['Class'] = 1
    ffp = pd.concat([globin_ffp, zincfinger_ffp], ignore_index=True)

    # Task 5 - classification
    print(classification(ffp))
