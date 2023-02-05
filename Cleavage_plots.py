#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  4 19:32:13 2023

@author: Amy
"""

from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import csv

'''
This function loads cleavages into Python

Input is a csv file with two columns. The first contains the Uniprot ID and the second contains the numerical position of P1' of the cleavage site
'''

def loadCleavages(csvfile):
    cleavagelist = []
    
    with open(csvfile, "r") as fileA:
        reader = csv.reader(fileA)
        next(reader)
        for row in reader:
            cleavagelist.append(row)
    return cleavagelist

'''
This function converts the list of cleavages to a dictionary where the keys are Uniprot IDs and the values are the numerical positions of P1'
'''
def cleavageDict(list):
    cleavage_dict = dict()
    
    for line in list:
        if line[0] in cleavage_dict:
            cleavage_dict[line[0]].append(line[1])
        else:
            cleavage_dict[line[0]] = [line[1]]
    return cleavage_dict

'''
The following lines of code load in a list of cleavages previously identified by subtiligase N terminomics and convert it to a dict
'''

previous = loadCleavages("Jurkat_etop_previous.csv")
SL_dict = cleavageDict(previous)

'''
The following lines of code load in a list of newly identified cleavages from CHOPPER and convert it to a dict
'''

infolist = loadCleavages('newcleavages.csv')
new_dict = cleavageDict(infolist)

'''
The following block finds the SwissProt record for all of the newly identified cleavages and the annotated domain boundaries, if any. It then finds the same information if there are cleavages that were previously identified in the same protein
'''

loclist=[]
cleavagelist = []
errorlist = []

f = open('SwissProt_human.xml', 'r')
for record in SeqIO.parse(f, 'uniprot-xml'):
    for key in new_dict.keys():
        if key in record.id:
            for ft in record.features:
                if "domain" in ft.qualifiers.get("type", ""):
                    domaintype = ft.qualifiers.get("description", "")
                    start = ft.location.nofuzzy_start + 1
                    end = ft.location.nofuzzy_end
                    length = len(record.seq)
                    name = record.id
                    loclist.append((name, domaintype, start, end, length))
                    cleavagelist = new_dict[key]
                    try:
                        previouslist = SL_dict[key]
                    except:
                        pass

#the next section plots cleavages and domain boundary information on a number line and writes a pdf of the plot for every protein in the list
                         
            try:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_xlim(0,10)
                ax.set_ylim(0,10)
                xmin = 0
                xmax = 10
                y = 8
                height = 0
                plt.hlines(y, -10, 20, color='w')
                plt.hlines(y, xmin, xmax)
                for i in range(len(loclist)):
                    start = int(loclist[i][2])
                    end = int(loclist[i][3])
                    length = int(loclist[i][4])
                    startval=(float(start)/float(length))*10
                    endval=(float(end)/float(length))*10
                    domainname=loclist[i][1]
                    middleval=np.mean((startval, endval))
                    plt.hlines(y=8, xmin=startval, xmax=endval, linewidth=16, color='b')
                    ax.annotate(domainname, xy=(startval,endval), xytext=(middleval, 7.3), rotation=90,fontsize=12, horizontalalignment='center', verticalalignment='top')        
                for j in range(len(cleavagelist)):
                    try:
                        arrowposition = (float(cleavagelist[j])/float(length))*10
                        plt.arrow(arrowposition, 9.2, 0.0, -0.4, fc="m", ec="m", head_width=0.15, head_length=0.3, linewidth=1)
                    except:
                        pass
                for k in range(len(previouslist)):
                    try:
                        arrowposition = (float(previouslist[k])/float(length))*10
                        plt.arrow(arrowposition, 9.2, 0.0, -0.4, fc="k", ec="k", head_width=0.15, head_length=0.3, linewidth=1)
                    except:
                        pass
                plt.text(xmin - 0.1, y, '1', horizontalalignment='right')
                plt.text(xmax + 0.1, y, length, horizontalalignment='left')
                plt.axis('off')
                outputfile = ''.join((loclist[i][0], '.pdf'))
                plt.savefig(outputfile, format="pdf")
            except:
                errorlist.append(key)
                #pass
            loclist = []
            cleavagelist = []
            previouslist = []
            