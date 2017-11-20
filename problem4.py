"""
Задача 4
Множественное выравнивание последовательностей
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


class LocalSequensAligner():
    def __init__(self, scoring_matrix, indel_penalty):
        self.a = {'A':0, 'T':1, 'G':2, 'C':3}
        self.sm = scoring_matrix
        self.indel = indel_penalty

    def align(self, seq1, seq2):
        alignment = np.zeros((len(seq1)+1,len(seq2)+1), dtype=np.int32)
        alignment_track = np.zeros((len(seq1)+1,len(seq2)+1), dtype=np.int32)
        max_alignment = 0
        max_alignment_pos = None
        for i in range(1, alignment.shape[0]):
            for j in range(1, alignment.shape[1]):
                options = np.array([alignment[i-1, j] + self.indel,
                    alignment[i, j-1] + self.indel,
                    alignment[i-1, j-1] + self.sm[self.a[seq1[i-1]], self.a[seq2[j-1]]]
                ])
                alignment_track[i,j] = np.argmax(options)
                alignment[i,j] = options[alignment_track[i,j]]

                if alignment[i,j] >= max_alignment:
                    max_alignment = alignment[i,j]
                    max_alignment_pos = [i,j]

        return alignment, alignment_track, max_alignment, max_alignment_pos


class GlobalSequenceAligner():
    pass


class Alignment():
    def __init__(self, note, seq1, seq2, aligner):
        self.name = note
        self.alignment, self.alignment_track, self.max_alignment, \
        self.max_alignment_pos = aligner.align(seq1, seq2)

        # Calculate alignment
        if self.max_alignment > 0:
            seq1_repr = []
            seq2_repr = []
            i, j = self.max_alignment_pos
            while i > 0 or j > 0:
                if self.alignment_track[i,j]==0:
                    seq1_repr.append(seq1[j-1])
                    seq2_repr.append('-')
                    i-=1
                if self.alignment_track[i,j]==1:
                    seq2_repr.append(seq2[j-1])
                    seq1_repr.append('-')
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
        else:
            self.seq1_repr = seq1
            self.seq2_repr = seq2

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


parser = argparse.ArgumentParser(description='Sorting of fasta file (dna).')
parser.add_argument('--file', type=str, help='path to .fasta file')
parser.add_argument('--scoring_matrix', type=str, help='path to similarity_matrix')
parser.add_argument('--indel_penalty', type=float, default=0.5, help='indel penalty')
parser.add_argument('--ref', type=str, help='reference string')
args = parser.parse_args()

scoring_matrix = read_array(args.scoring_matrix)
indel_penalty = args.indel_penalty
aligner = LocalSequensAligner(scoring_matrix, indel_penalty)

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

sorted_alignments = tree.in_order_traversal()
for a in reversed(sorted_alignments):
    print(a[1])
