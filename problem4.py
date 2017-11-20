"""
Задача 4
Множественное выравнивание последовательностей
"""

import sys, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import numpy as np


class FastaSeq:
    def __init__(self, type, name, seq):
        self.type = type
        self.name = name
        self.seq = seq
        if type=='dna':
            self.alphabet = alphabet_dna
        elif type=='protein':
            self.alphabet = alphabet_amino

    def __str__(self):
        return ">"+self.name+"\n"+self.seq+"\n"

    def __repr__(self):
        return ">"+self.name+"\n"+self.seq+"\n"


class LocalSequensAlignment():
    def __init__(scoring_matrix, indel_penalty):
        def.alph =




def read_array(file_path):
    with open(file_path) as f:
        content = f.readlines()
        arr = np.zeros((len(content), len(content[0].split(' '))))
        for i,line in enumerate(content):
            arr[i,] = [float(el) for el in line.split(' ')]
    return arr

print(read_array('./example_data/similarity_matrix.txt'))


parser = argparse.ArgumentParser(description='Sorting of fasta file (protein or dna).')
parser.add_argument('--file', type=str, help='path to .fasta file')
parser.add_argument('--score_file', type=str, help='path to similarity_matrix')
parser.add_argument('--indel_penalty', type=float, default=0.5 help='indel penalty')
args = parser.parse_args()

sequences = []
with open(args.file) as handle:
    for values in SimpleFastaParser(handle):
        sequences.append(FastaSeq(args.type, values[0], values[1]))

for seq in heapSort(sequences):
    print(seq)
