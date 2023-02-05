#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 19:20:16 2023

@author: Amy
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

'''
Read in N termini with cleaved biotin-clicked alkyne-2PCA modification. Starting from a comprehensive enhanced annotation Proteome Discoverer output peptide groups csv file, make a list of proteins and locations of the identified N terminus within the protein
for peptides with a user-defined N-terminal modification. 

Inputs required are filename, modification

Output is a two-column list containing Uniprot ID and numerical start site of the N-terminal peptide.
'''

def readModifiedNeoNTermini(file, modification): 
    filename = file
    mod = modification
    proteinandpos = []
   
    with open(filename, "r") as fileA:
        reader = csv.reader(fileA)
        next(reader)
        for row in reader:
            if mod in row[3]:
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

#count = countNonPrimePrimeAAs(nonprime)

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
Calculate z-scores given two input arrays (e.g., one that count amino acid occurence in modified peptides and one that counts amino acid
occurence in all peptides in a sample)

Required inputs are two arrays of the same size (sample and control). Output is an array of z-scores. This is also written to a file called
output.csv that will be overwritten everytime this function is used. Need to either rename the file or change the function if you want to keep
each file.
'''

def standardScore2(sample, control): #sample and control are arrays
    sampletotal = np.sum(sample, axis=0) #axis=0 sums the columns, axis=1 sums over rows
    samplefreq = sample/sampletotal
    controltotal = np.sum(control, axis=0)
    controlfreq = control/controltotal
    stderror = np.sqrt(controlfreq*(1-controlfreq)*((1/sampletotal)))
    zscore = ((samplefreq-controlfreq)/stderror)
    np.savetxt('output.csv', zscore, delimiter = ',', newline='\n') #saves a csv file with the calculated zscores
    return zscore


'''
Plot the z-scores generated by standardScore2 as a heatmap. Required inputs are an array of z-scores and a pdf filename. Output is a pdf
of the heatmap.
'''

def plotZscoreHeatmap(inputarray, pdfname):
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    maxvalue=np.nanmax(inputarray)
    minvalue=np.nanmin(inputarray)
    
    #set axis limits for heatmap so that 0 (i.e., no enrichment or deenrichment) is white
    if abs(maxvalue) > abs(minvalue):
        highlim=abs(maxvalue)
        lowlim=-highlim
    elif abs(minvalue) > abs(maxvalue):
        highlim=abs(minvalue)
        lowlim=-highlim
    
    #plot figure and save as pdf
    fig, ax = plt.subplots()
    ax = sns.heatmap(inputarray, vmin=lowlim, vmax=highlim, yticklabels=aminoacids, square=True, cmap='RdBu', clim=(-20,20))
    
    #rotate amino acid labels
    plt.yticks(rotation=0)
    
    #next two lines expand x and y limits so boxes aren't cut off at the top and bottom
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    
    #save heatmap as pdf
    fig.savefig(f"{pdfname}_fractionoccurence.pdf") #can change the name here


'''
trypsin

'''
E03_neoNterm = readModifiedNeoNTermini('E20220729-03.csv', "1xclicked DTT 2PCA [N-Term]")
E03_nonprime = findNonPrimeSeq(E03_neoNterm, 'SwissProt_human.xml')
E03_count = countNonPrimePrimeAAs(E03_nonprime)
E03_frequency = frequencyAAOccurence(E03_count)
E03_plotfrequency = plotFractionOccurence(E03_frequency, 'E20220729-03')


E05_neoNterm = readModifiedNeoNTermini('E20220729-05.csv', "1xclicked DTT 2PCA [N-Term]")
E05_nonprime = findNonPrimeSeq(E05_neoNterm, 'SwissProt_human.xml')
E05_count = countNonPrimePrimeAAs(E05_nonprime)
E05_frequency = frequencyAAOccurence(E05_count)
E05_plotfrequency = plotFractionOccurence(E05_frequency, 'E20220729-05')

E07_neoNterm = readModifiedNeoNTermini('E20220729-07.csv', "1xclicked DTT 2PCA [N-Term]")
E07_nonprime = findNonPrimeSeq(E07_neoNterm, 'SwissProt_human.xml')
E07_count = countNonPrimePrimeAAs(E07_nonprime)
E07_frequency = frequencyAAOccurence(E07_count)
E07_plotfrequency = plotFractionOccurence(E07_frequency, 'E20220729-07')

E09_neoNterm = readModifiedNeoNTermini('E20220729-09.csv', "1xclicked DTT 2PCA [N-Term]")
E09_nonprime = findNonPrimeSeq(E09_neoNterm, 'SwissProt_human.xml')
E09_count = countNonPrimePrimeAAs(E09_nonprime)
E09_frequency = frequencyAAOccurence(E09_count)
E09_plotfrequency = plotFractionOccurence(E09_frequency, 'E20220729-09')

total_etop_trypsin = E03_count + E05_count
total_DMSO_trypsin = E07_count + E09_count

zscores = standardScore2(total_etop_trypsin, total_DMSO_trypsin)
np.set_printoptions(precision=3, suppress=True)
zscores2 = np.nan_to_num(zscores, nan=0, posinf = 0)

heatmap = plotZscoreHeatmap(zscores2, 'etoposide_vs_DMSO')

#write frequencies to csv file to use in bar graphs, etc
   
np.savetxt("E03_frequency.csv", E03_frequency, delimiter = ",")
np.savetxt("E05_frequency.csv", E05_frequency, delimiter = ",")
np.savetxt("E07_frequency.csv", E07_frequency, delimiter = ",")
np.savetxt("E09_frequency.csv", E09_frequency, delimiter = ",")
