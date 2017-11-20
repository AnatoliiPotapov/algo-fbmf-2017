"""
Задача 5. Картирование ридов на референс при помощи полиномиальных хешей.

Картирование коротких прочтений
(англ. Short-Read Sequence Alignment, Short-Read Sequence Mapping) —
 биоинформатический метод анализа результатов секвенирования нового поколения,
 состоящий в определении позиций в референсном геноме или транскриптоме,
 откуда с наибольшей вероятностью могло быть получено каждое конкретное
 короткое прочтение. Обычно является первой стадией в обработке данных в случае,
 если известен геном исследуемого организма.
"""
import sys, argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import argparse

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

    def compare(self, seq):
        return compare(self.alphabet, self, seq)


class RollingHashBio:
    def __init__(self, reference, pattern_size, seq_type):
        if seq_type == 'dna':
            self.alphabet = alphabet_dna
            self.pow = 5
        elif seq_type == 'protein':
            self.alphabet = alphabet_amino
            self.pow = 23

        self.ref = reference
        self.ptrn = pattern_size
        self.hash = 0

        for i in range(0, pattern_size):
            self.hash += (self.alphabet[self.ref[i]]+1)*(self.pow**(self.ptrn - i -1))

        #start index of current window
        self.window_start = 0
        #end of index window
        self.window_end = self.ptrn

    def move_window(self):
        if self.window_end <= len(self.ref) - 1:
            #remove left letter from hash value
            self.hash -= (self.alphabet[self.ref[self.window_start]]+1)*self.pow**(self.ptrn-1)
            self.hash *= self.pow
            self.hash += self.alphabet[self.ref[self.window_end]]+1
            self.window_start += 1
            self.window_end += 1

    def window_text(self):
        return self.ref[self.window_start:self.window_end]

    def digest(self):
        return self.hash


def single_pattern_rabin_karp(ref, reed, seq_type):
    if ref== "" or reed == "":
        return None
    if len(reed) > len(ref):
        return None

    all_occurrences = []

    ref_rolling_hash = RollingHashBio(ref, len(reed), seq_type)
    reed_hash = RollingHashBio(reed, len(reed), seq_type)

    for i in range(len(ref) - len(reed) + 1):
        if ref_rolling_hash.hash == reed_hash.hash:
            if ref_rolling_hash.window_text() == reed:
                all_occurrences.append(i)
        ref_rolling_hash.move_window()

    if len(all_occurrences)>0:
        return all_occurrences
    return None


def multiple_pattern_rabin_karp(ref, reeds, seq_type):
    pass


parser = argparse.ArgumentParser(description='Finding all occurrences of reeds \
from .fasta file in reference sequence (protein or dna).')
parser.add_argument('--type', type=str, help='type of sequences')
parser.add_argument('--file', type=str, help='path to .fasta file')
parser.add_argument('--ref', type=str, help='reference sequence')
args = parser.parse_args()

sequences = []
with open(args.file) as handle:
    for values in SimpleFastaParser(handle):
        sequences.append(FastaSeq(args.type, values[0], values[1]))

test_reference = args.ref
reed = sequences[0].seq
print(single_pattern_rabin_karp(test_reference, reed, 'dna'))
