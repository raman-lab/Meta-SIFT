#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison
# Feb 2020 (updated Dec 2021)


# Usage: $ motif_finder_tool.py -i <input_proteins> -m <matrix_table> [options]

import argparse
import os
import subprocess
import numpy as np
import motif_finder_tool_modules as cy


# Publication setttings
# -i <one of the protein datasets>
# -m <f-score table>

# For 6-mer analysis
# -c 1
# -n 6
# -e 1e-50
# -t 0.8

# For 10-mer analysis
# -c 1
# -n 10
# -e 1e-5
# -t 0.045

motif_finder_path = str(os.path.dirname(os.path.abspath(__file__)))
motif_finder = argparse.ArgumentParser(description='Identifies proteins with similar motifs to the input matrix table.')
motif_finder.add_argument('-i', type=str, nargs=1, required=True, help='input proteins fasta file.')
motif_finder.add_argument('-m', type=str, nargs=1, required=True, help='input maxtrix table of raw F-score values. No header. Column A = WT amino acid letters. Columns B-U = F-scores of each of the 20 possible amino acids (alphabetical order). No stop codon scores.')
#
motif_finder.add_argument('-o', type=str, nargs=1, default='not set', help='output results table of scores [default=input.motif-scores.tsv].')
motif_finder.add_argument('-a', type=str, nargs=1, default='not set', help='output results table of motif hits summary [default=input.motif-summary.tsv].')
motif_finder.add_argument('-k', type=str, nargs=1, default='not set', help='output results table of motif sequences [default=input.motif-seqs.tsv].')
motif_finder.add_argument('-n', type=str, nargs=1, default='6', help='nmer size.')
motif_finder.add_argument('-c', type=str, nargs=1, default='not set', help='arbitrary score cutoff per nmer hit.')
motif_finder.add_argument('-e', type=float, nargs=1, default='1e-50', help='maximum threshold of evalue to report.')
motif_finder.add_argument('-s', type=str, nargs=1, default='0', help='minimum protein size. Must be >= length of matrix Column A.')
motif_finder.add_argument('-l', type=str, nargs=1, default='2000', help='maximum protein size.')
motif_finder.add_argument('-t', type=str, nargs=1, default='cutoff value', help='minimum threshold of normalized scores to report.')

#
args = motif_finder.parse_args()
#
try:
    matrix = str(args.m[0]).rsplit("/",1)[1]
except Exception:
    matrix = str(args.m[0])
strain = str(matrix).rsplit(".",1)[0]
#
if type(args.c) is list:
    cutoff = float(args.c[0])
else:
    cutoff = 1
#
if type(args.n) is list:
    length = int(args.n[0])
else:
    length = int(args.n)
#
if type(args.o) is list:
    output = str(args.o[0])
else:
    output = str(args.i[0]).rsplit(".",1)[0]+'.' + str(strain) + '.c-' + str(cutoff) + '_n-' + str(length) + '.motif-scores.tsv'
#
if type(args.a) is list:
    outfile = str(args.a[0])
else:
    outfile = str(args.i[0]).rsplit(".",1)[0]+'.' + str(strain) + '.c-' + str(cutoff) + '_n-' + str(length) + '.motif-summary.tsv'
#
if type(args.k) is list:
    outseqs = str(args.k[0])
else:
    outseqs = str(args.i[0]).rsplit(".",1)[0]+'.' + str(strain) + '.c-' + str(cutoff) + '_n-' + str(length) + '.motif-seqs.tsv'
#
#
if type(args.s) is list:
    min_size = int(args.s[0])
else:
    min_size = int(args.s)
#
if type(args.l) is list:
    max_size = int(args.l[0])
else:
    max_size = int(args.l)
#
if type(args.t) is list:
    threshold = float(args.t[0])
else:
    threshold = 0.8

if type(args.e) is list:
    evalue = float(args.e[0])
else:
    evalue = 1e-50

aa_table = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] # used for all-strains and legacy tables

def fasta_parse(infile):
    with open(infile, 'r') as fasta:
        for line in fasta:
            line = line.strip("\n")
            if line[0] == '>':
                try:
                    yield header, ''.join(seq)
                except NameError:
                    pass # first line
                seq = []
                header = line[1:]
            else:
                seq.append(line)

        # last one
        yield header, ''.join(seq)

############################################################### read in table and generate final "transform" list
with open(str(args.m[0]), 'r') as mot:
    mot = mot.read().split("\n")
    len_mot = len(mot[0].split("\t"))

if min_size <= len(mot):
    min_size = len(mot)

w = 0
transform = [] # columns of each amino acid in the motif that have been transformed
                # transformation is sum of column and then each item divided by the total sum of the column
                # this normalizes each column to a total value of 1
t_count = 0
k = 0 # used with mot_dict and mot_list
mot_dict = {} # used to count motif hits (store hit counts)
mot_list = [] # used to count motif hits (store amino acids)

while w <= len(mot)-length:
    motif = []
    wt = ''
    for item in mot[w:w+length]:
        mot_mer = item.split("\t")
        motif.append(mot_mer[1:])
        wt += str(mot_mer[0])
    t_count += 1
    mot_list.append(str(wt))
    mot_dict.update({k:0})
    k += 1

    t = 0
    temp = []
    for item in motif:
        for nums in item:
            temp.append(float(nums)) #*** this change was made to eliminate normalization per column
        t += 1
    transform.append(temp) # transform is the transformed (normalized) matrix scores

    w += 1 # next window

zip_table = []
for item in transform:
    zip_table.append(dict(zip(aa_table,item))) # zip table contains the score for each of the 20 amino acids for each of the nmer amino acids in the window

zip_table = [zip_table[i:i+length] for i in range(0,len(zip_table)-length+1)]

len_trans = len(transform)
################################################## start the scan
with open(str(output), 'w') as summary, open(str(outfile), 'w') as count_file, open(str(outseqs), 'w') as outseqs:
    summary.write("protein\tnmer value\tnmer count\tseq length\twindows scored\twindows total\tunique motifs\ttotal score\tnormalized score\tevalue\n")
    outseqs.write("protein\tmotif table\tscore\tstart\tlength\tWT_motif\tquery_motif\n")
    count_file.write("WT_motif")
    for name in mot_list:
        count_file.write("\t" + str(name))
    count_file.write("\n")
    t = 0
    total = 0
    len_trans_length = len_trans-length
    for name,seq_aa in fasta_parse(args.i[0]):
        seq_aa = seq_aa.replace('\n','').strip('*').upper()
        if ' # ' in name:
            name = name.split(' # ')[0]
        elif '\t' in name:
            name = name.split('\t')[0]
        name = name.replace(' ','~~')
        Xs = seq_aa.count('X')
        seqs_list = []
        length_check,l_aa,l_aa_Xs = cy.min_max_size(seq_aa, Xs, min_size, max_size)

        if length_check:
            score_max = 0
            windows = 0
            win_check = 0
            zeros = 0
            kmers = [seq_aa[i:i+length] for i in range(0,l_aa-length+1)]
            kmers = [k for k in kmers if 'X' not in k]
            for i,nmer in enumerate(kmers):
                master = np.zeros(len_trans)
                # iterate through all the matrices
                for t,mot_table in enumerate(zip_table):
                    
                    aa_score = 0
                    score_check = []
                    for idx,aa in enumerate(nmer):
                        table = mot_table[idx]
                        score = table[aa]
                        score_check.append(score) # used in "if" statement below with cutoff

                    l_sc = len(score_check)
                    len_check = cy.score_check_len(l_sc, length)
                    if len_check:
                        win_check = cy.win_check_add(win_check) # this is how many windows were checked
                        score_check.sort()
                        filter_scores_check,aa_score = cy.filter_scores(score_check, length, cutoff)
                        if filter_scores_check:
                            name_replace = name.replace('~~',' ')
                            master[t] = aa_score
                            score_check = []
                            windows = cy.windows_add(windows) # this is how many windows were scored
                            mot_dict[t] += 1 # counting the motif that it hit
                            seqs_list.append(f'{name_replace}\t{t}\t{aa_score}\t{i}\t{l_aa}\t{mot_list[t]}\t{nmer}\n')
                        else:
                            score_check = []
                    else:
                        score_check = []
                
                m = np.ndarray.max(master)
                score_max = cy.get_score_max(score_max, m, length)

            md_values = list(mot_dict.values())
            zeros = cy.get_zeros(md_values)
            final, e_val = cy.get_final_scores(score_max, l_aa_Xs, windows, win_check)

            check_e_f = cy.evalue_final(windows, final, threshold, e_val, evalue)
            if check_e_f:

                for seq in seqs_list:
                    outseqs.write(seq)
                summary.write(f'{name_replace}\t{length}\t{t_count}\t{l_aa}\t{windows}\t{win_check}\t{zeros}\t{score_max}\t{final}\t{e_val}\n')

                count_file.write(name_replace)
                sum_p = sum(md_values)
                for p in md_values:
                    new = cy.mot_dict_math(p, sum_p)
                    count_file.write(f"\t{new}")
                count_file.write("\n")
            for key in mot_dict.keys():
                mot_dict[key] = 0 # reset all the dictionary keys back to zero for the next protein


subprocess.run("mkdir " + str(args.i[0]).rsplit(".",1)[0] + '_' + str(strain) + '_c-' + str(cutoff) + '_n-' + str(length) + "_motif-finder", shell=True)
subprocess.run("mv " + str(args.i[0]).rsplit(".",1)[0] + '.' + str(strain) + '.c-' + str(cutoff) + '_n-' + str(length) + ".* " + str(args.i[0]).rsplit(".",1)[0] + '_' + str(strain) + '_c-' + str(cutoff) + '_n-' + str(length) + "_motif-finder", shell=True)


#
#
#
