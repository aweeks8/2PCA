#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 22:38:21 2022

@author: Amy
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import re

'''
Starting from a comprehensive enhanced annotation Proteome Discoverer output peptide groups csv file, make a list of proteins and locations of the identified N terminus within the protein
for peptides with a user-defined N-terminal modification. Ignores N termini with same specificity as the digest protease.

Inputs required are filename, modification (e.g. '2PCA' for PICS2, 'Thioacylation' for PICS), and identity of digest protease (chymo, gluc, or tryp). 

Output is a two-column list containing Uniprot ID and numerical start site of the N-terminal peptide.
'''

def readModifiedNeoNTermini(file, modification, libraryprotease): #ignores N termini with same specificity as library generation protease
    filename = file
    mod = modification
    proteinandpos = []
    
    if libraryprotease == 'chymo':
        prospecificity = "[FYWML]"
    if libraryprotease == 'gluc':
        prospecificity = "[DE]"
    if libraryprotease == 'tryp':
        prospecificity = "[KR]"
    if libraryprotease == 'chymo_low':
        prospecificity = "[FYWMLH]"
    
    with open(filename, "r") as fileA:
        reader = csv.reader(fileA)
        next(reader)
        for row in reader:
            if mod in row[3]:
                P1 = row[2].split(".")[0][1]
                match = re.search(prospecificity, P1)
                if match:
                    pass
                else:
                    proteinandpos.append(row[10])
    return proteinandpos

'''
This function finds the prime and nonprime sequences given a protein ID and position of residue on the C-terminal side of scissile bond.

Inputs required are 1)a list of Uniprot IDs and positions (output of readModifiedNeoNTermini) and the name of a Uniprot XML file.

Output is a list of sequences with four residues on the prime side and four residues on the nonprime side with the cleavage site indicated by a '|'.
'''

def findNonPrimeSeq(modifiedproteins, uniprot_file):

    infolist = []
    
    for row in modifiedproteins:
        proteinnameend = row.index(" ")
        print(row)
        protein = row[0:proteinnameend]
        peptidestartbeg = row.index("[")
        peptidestartend = row.index("-")
        peptidestart = int(row[peptidestartbeg+1:peptidestartend])
        peptideendbeg = row.index("-")
        peptideendend =row.index("]")
        peptideend = int(row[peptideendbeg+1:peptideendend])
        infolist.append((protein, peptidestart, peptideend))
        
    
    #use the proteinandposlist to find the position of the observed peptide within protein sequences in the Uniprot xml file and to pull out the preceding amino acids
    peptidelist = []
    f = open(uniprot_file, 'r')
    for record in SeqIO.parse(f, 'uniprot-xml'):
        for i in range(len(infolist)):
    #the first block applies to only peptides that begin at position 5 or later (because this will give us a valid position in the sequence)
            if infolist[i][1] >= 5:
                if infolist[i][0] in record.id:
                    print(infolist[i][0], record.seq[infolist[i][1]-1:infolist[i][2]])
                    print(record.seq[infolist[i][1]-5:infolist[i][1]-1])
                    print(record.seq[infolist[i][1]-5:infolist[i][1]-1]+"|"+record.seq[infolist[i][1]-1:infolist[i][2]])
                    peptidelist.append(record.seq[infolist[i][1]-5:infolist[i][1]-1]+"|"+record.seq[infolist[i][1]-1:infolist[i][2]])
    #the else block applies to peptides with a start position <5, so we can get as many amino acids as there are before the cleavage site at the beginning of the protein                
            else:
                if infolist[i][0] in record.id:
                    print(infolist[i][0], record.seq[infolist[i][1]-1:infolist[i][2]])
                    print(record.seq[0:infolist[i][1]-1])
                    print(record.seq[0:infolist[i][1]-1]+"|"+record.seq[infolist[i][1]-1:infolist[i][2]])
                    peptidelist.append(record.seq[0:infolist[i][1]-1]+"|"+record.seq[infolist[i][1]-1:infolist[i][2]])
                    
    return peptidelist

'''
This function accepts the output of findNonPrimeSeq as input and counts the number of each amino acid at each position

Output is a 8 x 20 numpy array with amino acids counts. The 8 columns correspond to positions P4-P4' and the 20 rows correspond to the 20 amino acids in alphabetical order.
'''

def countNonPrimePrimeAAs(peptidelist):

    #Set up dictionary that encodes amino acids as row numbers
                    
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
    val = 0
    aminoDict=dict()
        
    for AA in aminoacids:
        aminoDict[AA] = val
        val += 1
    
    #make a 20 x 9 array populated with zeroes where we will eventually plot the amino acid counts in each position
    peptidecount = np.zeros((20,9)) 
    
    #count the amino acids in each position on the nonprime (before the "|") and prime (after the "|") sides
    for row in peptidelist:
        try:
            pepseq=str(row) #I'm not sure if this would be necessary with the newest version of Biopython
            cleavagemarker = pepseq.index("|")
            for i in range((cleavagemarker-4), (cleavagemarker+5)):
                try:
                    peptidecount[aminoDict[row[i]]][i-(cleavagemarker-4)] +=1
                except:
                    pass
        except:
            pass
    
    #delete the row of the peptide count array that corresponds to "|", the cleavage marker
    cleavagemarkerremoved = np.delete(peptidecount, 4, 1)

    return cleavagemarkerremoved

'''
This function accepts the output of countNonPrimePrimeAAs and converts it to the frequency of each amino acid in each position instead of count. Does this by dividing by the number of observations in each column

Output is an 8 x 20 array of frequencies. The 8 columns correspond to positions P4-P4' and the 20 rows correspond to the 20 amino acids in alphabetical order.
'''

def frequencyAAOccurence(inputarray):
    fractionoccurence = 100*(inputarray/np.sum(inputarray, axis = 0))
    return fractionoccurence

'''
This function plots the frequency array as a heatmap. Required inputs are the name of the array and a title for saving an image of the array as a PDF
'''

def plotFractionOccurence(inputarray, title):
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    xticklabelslist = ['-4', '-3', '-2', '-1', '0', '1', '2', '3', '4']
    fig, ax = plt.subplots()
    sns.color_palette(palette="viridis")
    ax = sns.heatmap(inputarray, vmin=0, vmax=100, yticklabels=aminoacids, xticklabels=xticklabelslist, square=True, cmap="viridis")
    cbar = ax.collections[0].colorbar
    cbar.set_label('amino acid occurence (%)', rotation=270, labelpad=20)
    
        #rotate y-axis labels
    plt.yticks(rotation=0)
    
    #next two lines expand x and y limits so boxes aren't cut off at the top and bottom
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    
    #move x-axis labels to top and rotate them
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)
    
    #save figure as a pdf file
    fig.savefig(f"{title}_fractionoccurence.pdf")

    return
'''
The lines below show the use of the above functions to analyze Kex2 PICS2 data.
'''
#analysis
    
Kex2_neo_termini_chymo = readModifiedNeoNTermini('E20221129-02.csv', "1xreduced biotinSS2PCA [N-Term]", 'chymo')

nonprime_Kex2_chymo = findNonPrimeSeq(Kex2_neo_termini_chymo, 'SwissProt_Ecoli_K12.xml')

Kex2_count_chymo = countNonPrimePrimeAAs(nonprime_Kex2_chymo)

Kex2_fractionoccurence_chymo = frequencyAAOccurence(Kex2_count_chymo)

Kex2_plot_chymo = plotFractionOccurence(Kex2_fractionoccurence_chymo, 'Kex2_chymo_Ecoli_fraction_PICS2')

#Write N-terminal sequences to csv
with open('Kex2_nonprime_prime_chymo.csv', "w") as fileC:
	writer = csv.writer(fileC, lineterminator = '\n')
	writer.writerows(nonprime_Kex2_chymo )
    
    
Kex2_neo_termini_gluc = readModifiedNeoNTermini('E20221129-04.csv', "1xreduced biotinSS2PCA [N-Term]", 'gluc')

nonprime_Kex2_gluc = findNonPrimeSeq(Kex2_neo_termini_gluc, 'SwissProt_Ecoli_K12.xml')

Kex2_count_gluc = countNonPrimePrimeAAs(nonprime_Kex2_gluc)

Kex2_fractionoccurence_gluc = frequencyAAOccurence(Kex2_count_gluc)

Kex2_plot_gluc = plotFractionOccurence(Kex2_fractionoccurence_gluc, 'Kex2_gluc_Ecoli_fraction_PICS2')

#Write N-terminal sequences to csv
with open('Kex2_nonprime_prime_gluc.csv', "w") as fileC:
	writer = csv.writer(fileC, lineterminator = '\n')
	writer.writerows(nonprime_Kex2_gluc)

#Combine chymotrypsin and GluC datasets 
    
combined = nonprime_Kex2_chymo + nonprime_Kex2_gluc

combined_count = countNonPrimePrimeAAs(combined)
Kex2_fractionoccurence_combined = frequencyAAOccurence(combined_count)
Kex2_plot_combined = plotFractionOccurence(Kex2_fractionoccurence_chymo, 'Kex2_combined_Ecoli_fraction_PICS2')

with open('Kex2_nonprime_prime_combined.csv', "w") as fileC:
	writer = csv.writer(fileC, lineterminator = '\n')
	writer.writerows(combined)

