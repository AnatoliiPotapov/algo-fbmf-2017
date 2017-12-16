"""
Задача 4
Множественное выравнивание последовательностей

На вход программе подаётся референсная строка и fasta-файл, в котором содержатся риды.
Необходимо отсортировать риды по величине оценки выравнивания с помощью двоичного дерева.

Решение: для всех ридов вычисляется alignment_score, при помощи локального(глобального)
выравнивания, далее риды сортируются по величине alignment_score при помощи двоичного
дерева.
"""

import sys, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
import numpy as np
from binary_tree import Tree, Node

alphabet_dna = {'_':-1,'A':0,'C':1,'G':2,'T':3}
alphabet_amino = {'_':-1,'A':0,'R':1,'N':2,'D':3,'B':4,'C':5,'E':6,'Q':7,
    'Z':8,'G':9,'H':10,'I':11,'L':12,'K':13,'M':14,'F':15,'P':16,'S':17,
    'T':18,'W':19,'Y':20,'V':21}

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

# Works only for dna!
class SequenseAligner():
    def __init__(self, scoring_matrix, indel_penalty, mode = 'global'):
        self.a = {'A':0, 'T':1, 'G':2, 'C':3}
        self.sm = scoring_matrix
        self.indel = indel_penalty
        self.mode = mode

    def align(self, seq1, seq2):
        alignment = np.zeros((len(seq1)+1,len(seq2)+1), dtype=np.float32)

        if self.mode == 'global':
            alignment[:,0] = np.array([float(self.indel) * i for i in range(len(seq1)+1)])
            alignment[0,:] = np.array([float(self.indel) * i for i in range(len(seq2)+1)])

        alignment_track = np.zeros((len(seq1)+1,len(seq2)+1), dtype=np.int32)
        max_alignment = 0
        max_alignment_pos = None
        for i in range(1, alignment.shape[0]):
            for j in range(1, alignment.shape[1]):
                options = np.array([
                    alignment[i-1, j] + self.indel,
                    alignment[i, j-1] + self.indel,
                    alignment[i-1, j-1] + self.sm[self.a[seq1[i-1]], self.a[seq2[j-1]]]
                ])

                if self.mode == 'local':
                    options = np.append(options, 0)

                alignment_track[i,j] = np.argmax(options)
                alignment[i,j] = options[alignment_track[i,j]]

                if alignment[i,j] >= max_alignment:
                    max_alignment = alignment[i,j]
                    max_alignment_pos = [i,j]

        if self.mode == 'global':
            max_alignment = alignment[-1,-1]
            max_alignment_pos = [len(seq1), len(seq2)]

        return alignment, alignment_track, max_alignment, max_alignment_pos

# Testing global alignment
#scoring_matrix = np.array([[1,-1,-1,-1], [-1,1,-1,-1], [-1,-1,1,-1], [-1,-1,-1,1]])
#aligner = SequenseAligner(scoring_matrix, -2, mode='global')
#true_answ = np.array([[ 0, -2, -4, -6], [-2, -1, -3, -5], [-4, -1, -2, -4], [-6, -3,  0, -2], [-8, -5, -2,  1]])
#print(aligner.align('ATCG','TCG')[0] ==  true_answ)

class Alignment():
    def __init__(self, note, seq1, seq2, aligner):
        self.name = note
        self.alignment, self.alignment_track, self.max_alignment, \
        self.max_alignment_pos = aligner.align(seq1, seq2)

        # Calculate alignment
        seq1_repr = []
        seq2_repr = []
        i, j = self.max_alignment_pos

        l = len(seq1)
        if aligner.mode == 'local':
            while l > i:
                l-=1
                seq1_repr.append(seq1[l])
                seq2_repr.append(' ')

        while i > 0 or j > 0:
            indel_symb = '-'
            if aligner.mode == 'local':
                if i==0 or j==0 or i==len(seq1) or j==len(seq2):
                    indel_symb = ' '
                else:
                    indel_symb = '-'

            if self.alignment_track[i,j]==0:
                seq1_repr.append(seq1[i-1])
                seq2_repr.append(indel_symb)
                i-=1
            if self.alignment_track[i,j]==1:
                seq2_repr.append(seq2[j-1])
                seq1_repr.append(indel_symb)
                j-=1
            if self.alignment_track[i,j]==2:
                seq1_repr.append(seq1[i-1])
                if seq1[i-1]==seq2[j-1]:
                    seq2_repr.append(seq2[j-1])
                else:
                    seq2_repr.append(seq2[j-1].lower())
                i-=1
                j-=1

        self.seq1_repr = ''.join(reversed(seq1_repr))
        self.seq2_repr = ''.join(reversed(seq2_repr))


    def get_score(self):
        return self.max_alignment

    def __str__(self):
        return self.name+'\n'+'score:'+str(self.max_alignment)+'\n'+\
            self.seq1_repr+'\n'+self.seq2_repr+'\n'

    def __repr__(self):
        return self.name+'\n'+'score:'+str(self.max_alignment)+'\n'+\
            self.seq1_repr+'\n'+self.seq2_repr+'\n'

def read_array(file_path):
    with open(file_path) as f:
        content = f.readlines()
        arr = np.zeros((len(content), len(content[0].split(' '))))
        for i,line in enumerate(content):
            arr[i,] = [float(el) for el in line.split(' ')]
    return arr

#scoring_matrix = read_array('./example_data/similarity_matrix.txt')
#aligner = SequenseAligner(scoring_matrix, -0.5, mode='global')
#aligner.align('ATTTCGCACTTGG', 'GGGG')

parser = argparse.ArgumentParser(description='Sorting of fasta file (dna).')
parser.add_argument('--file', type=str, help='path to .fasta file')
parser.add_argument('--scoring_matrix', type=str, help='path to similarity_matrix')
parser.add_argument('--indel_penalty', type=float, default=-0.5, help='indel penalty')
parser.add_argument('--mode', type=str, default='global', help='type of alignment')
parser.add_argument('--ref', type=str, help='reference string')
args = parser.parse_args()

scoring_matrix = read_array(args.scoring_matrix)
indel_penalty = args.indel_penalty
aligner = SequenseAligner(scoring_matrix, indel_penalty, mode=args.mode)

sequences = []
with open(args.file) as handle:
    for values in SimpleFastaParser(handle):
        sequences.append(FastaSeq('dna', values[0], values[1]))

alignments = []
for seq in sequences:
    alignments.append(Alignment(seq.name, args.ref, seq.seq, aligner))

tree = Tree(need_payload=True)
for a in alignments:
    tree.add(Node(a.get_score(), a))

print('\n')

sorted_alignments = tree.in_order_traversal()
for a in reversed(sorted_alignments):
    print(a[1])
