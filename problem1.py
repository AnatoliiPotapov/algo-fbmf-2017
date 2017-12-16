"""
Задача 1

Сортировка входящего fasta-файла по лексикографическому порядку с помощью кучи.

Решение:
Читаем риды и сортируем кучей, собственно :)
"""

import sys, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse
alphabet_dna = {'_':-1,'A':0,'C':1,'G':2,'T':3}
alphabet_amino = {'_':-1,'A':0,'R':1,'N':2,'D':3,'B':4,'C':5,'E':6,'Q':7,
    'Z':8,'G':9,'H':10,'I':11,'L':12,'K':13,'M':14,'F':15,'P':16,'S':17,
    'T':18,'W':19,'Y':20,'V':21}

def compare(alphabet, str0, str1):
    cmp0=str0.seq+'_'
    cmp1=str1.seq+'_'
    ind=0
    while True:
        if cmp0[ind]>cmp1[ind]:
            return 0
        elif cmp0[ind]<cmp1[ind]:
            return 1
        elif cmp0[ind]=='_' and cmp1[ind]=='_':
            return 0
        ind+=1

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

    def compare(self, seq):
        return compare(self.alphabet, self, seq)

def heapSort(sequence):
    def swap_items(index1, index2):
        if sequence[index1].compare(sequence[index2])==1:                               # !
            sequence[index1], sequence[index2] = sequence[index2], sequence[index1]

    def sift_down(parent, limit):
        while True:
            child = (parent + 1) << 1 # То же, что и parent * 2 + 2
            if child < limit:
                if sequence[child].compare(sequence[child - 1])==1:                       # !
                    child -= 1
                swap_items(parent, child)
                parent = child
            else:
                break

    # Тело функции heap_sort
    length = len(sequence)
    # Формирование первичной пирамиды
    for index in range((length >> 1) - 1, -1, -1):
        sift_down(index, length)
    # Окончательное упорядочение
    for index in range(length - 1, 0, -1):
        swap_items(index, 0)
        sift_down(0, index)

    return sequence

parser = argparse.ArgumentParser(description='Sorting of fasta file (protein or dna).')
parser.add_argument('--type', type=str, help='type of sequences')
parser.add_argument('--file', type=str, help='path to .fasta file')
args = parser.parse_args()

sequences = []
with open(args.file) as handle:
    for values in SimpleFastaParser(handle):
        sequences.append(FastaSeq(args.type, values[0], values[1]))

for seq in heapSort(sequences):
    print(seq)
