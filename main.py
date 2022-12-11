#!/usr/bin/env python

# coding: utf-8
#  Copyright (c) 2022. Joshua Khorsandi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
import seaborn as sns
import subprocess
import math
from Bio.Align import AlignInfo
import collections

# loads the Acinetobacter, Chlamydia data into the pandas df
df_a = pd.read_csv('/Users/joshuakhorsandi/Desktop/IntroComp/Final/Data/AcinetobacterDNARegulator.tsv', sep='\t')

# Writes the converted sequenced back to an output fasta file
# converts sequences into SeqRecord Objects
sequences = []
for index, row in df_a.iterrows():
    sequences.append(
        SeqRecord(Seq(row['Sequence']), id=df_a.at[index, 'Entry'], description=df_a.at[index, 'Protein names'],
                  name="DNA-Binding Regulator"))

# pad the sequences in the list
max_len = max([len(s.seq) for s in sequences])
GAPs = "-"
for seq in sequences:
    padding = GAPs * (max_len - len(seq.seq))  # creating the padding string
    seq.seq += padding

# Write the sequence objects to a fasta file
with open("/Users/joshuakhorsandi/Desktop/IntroComp/Final/Data/Out/out.fasta", 'w') as file:
    SeqIO.write(sequences, file, "fasta")

# Uses clust to get the MSA from command line
# command = "python Util/clustalo.py " \
#           "--email jkhorsa1@jh.edu " \
#           "--stype protein --outfile alignment " \
#           "--sequence Data/Out/out.fasta"
# subprocess.run(command, shell=True)

# Get the aligned sequences
with open(
        '/Users/joshuakhorsandi/Desktop/IntroComp/Final/Data/Out/alignment.clustal_num.aln-clustal_num.clustal_num') as a_file:
    msa = AlignIO.read(a_file, 'clustal')


# ------ ANALYSIS ------- #
def shannon_entropy(sequence):
    """
    Calculates the shannon entropy
    @param sequence: iterable sequence object
    @return: shannon entropy for the sequence
    """
    h = 0
    m = len(sequence)
    codes = collections.Counter([tmp_base for tmp_base in sequence])
    for code in codes:
        n_i = codes[code]
        p_i = n_i / m
        entropy = p_i * math.log(p_i, 2)
        h += entropy
    return (-1) * h


# gets the entropy dataframe for plotting
P = {'ID': [],
     'Sequence': [],
     'Entropy': []}
s_info = AlignInfo.SummaryInfo(msa)
c_seq = s_info.gap_consensus()
for i in range(len(msa)):
    P['ID'].append(msa[i].id)
    P['Sequence'].append(msa[i].seq)
    P['Entropy'].append(shannon_entropy(msa[i].seq))
df_entropy = pd.DataFrame(P)
df_entropy = df_entropy[~(df_entropy == 0).any(axis=1)]
df_sorted = df_entropy.sort_values('Entropy')

# plots the data
sns.set_theme(style="whitegrid")
plt.figure(num=None, figsize=(20, 18), dpi=80, facecolor='w', edgecolor='r')
sns.barplot(y='ID', x="Entropy", data=df_sorted.head(40), orient='h')
sns.despine(left=True, bottom=True)
plt.title("Entropy vs Protein ID", fontsize=16)
plt.ylabel('Protein ID', fontsize=16)
plt.xlabel('Shannon Entropy', fontsize=16)
plt.show()

df_sorted = df_sorted.reset_index(drop=True)
plt.scatter(x=df_sorted.index, y=df_sorted['Entropy'], c='b')
plt.xlabel("Index")
plt.ylabel("Shannon Entropy")
plt.title("Shannon Entropy of All Proteins")
plt.show()
