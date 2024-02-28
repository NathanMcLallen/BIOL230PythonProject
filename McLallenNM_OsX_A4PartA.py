#!/usr/bin/env python
#Python program that runs through the files in 'YeastGenes' folder and finds the ratio of acidic to basic codons
# finds the number of trinucleotide repeats
# Same dependencies as the original myGeneParser

import fileinput
import os
import os.path
from os import path
import glob

#Reads in seqeunces as elements in list
# header holds the filename in position corresponding to position in sequence
header = list()
sequence = list()

# main program
def main():
    tot_file = 0
    path_foldername = "C:\\Users\\gabsbi\\Desktop\\code-examples\\other\\YeastGenes"
    foldername = 'YeastGenes'
    for filename in os.listdir(foldername):
        #print(filename)
        tot_file += 1
        temp = ""
        my_path = path.join(foldername, filename)
        firstLine = True
        for line in fileinput.input(files = (my_path)):
            if firstLine != True:
                temp += line
            firstLine = False
        filename = filename[:-4]
        header.append(filename)
        sequence.append(temp)
        #myFile = open("C:\\Users\\gabsbi\\Desktop\\code-examples\\other\\Yeast_RNAseq\\Nagalakshmi_2008_5UTRs_V64.gff3")
        #print(myFile.name)
        #print(myFile.readlines())
    #average = 0

    for i in range(len(sequence)):
       # average += gcContent(sequence[i], i)
       mRNA = transcribe(sequence[i])
       print(header[i])
      # print(mRNA) 
       acidOrBase(mRNA)
       findTriRepeats(sequence[i])
       print ("\n \n")
        #myFile = open("C:\\Users\\gabsbi\\Desktop\\code-examples\\other\\Yeast_RNAseq\\Nagalakshmi_2008_5UTRs_V64.gff3")
        #print(myFile.name)
        #print(myFile.readlines())
   # average /= tot_file
   # print("\nAverage GC Content is: ", round(average, 1))
    #myFile = open("C:\\Users\\gabsbi\\Desktop\\code-examples\\other\\Yeast_RNAseq\\Nagalakshmi_2008_5UTRs_V64.gff3")
    #print(myFile.name)
    #print(myFile.readlines())


# Transcribe sequence to mRNA. I didn't remove this function so that I could edit it and use the transcribed codons for one of my functions, I hope that's OK.
def transcribe(seq):
   # print("\nTranscribing.. \n")
   # seq = seq[::-1]
    for i in range(len(seq)):
        if seq[i] == "A":
            seq = seq[:i] + "U" + seq[i+1:]
        elif seq[i] == "T":
            seq = seq[:i] + "A" + seq[i+1:]
        elif seq[i] == "G":
            seq = seq[:i] + "C" + seq[i+1:]
        elif seq[i] == "C":
            seq = seq[:i] + "G" + seq[i+1:]
    return seq
        

# Counts number of acidic and basic codons and prints ratio
def acidOrBase(RNASeq):
    acidCount = 0
    baseCount = 0
    for i in range (0, len(RNASeq), 3):
        #print (RNASeq[i])
        codon = RNASeq[i:i + 3]
        #print (codon)
        firstTwo = RNASeq[i:i + 2]
        if (firstTwo == "GA"):
            acidCount += 1
        elif (firstTwo == "CG" or codon == "AGA" or codon == "AGG" or codon == "AAA" or codon == "AAG" or codon == "CAU" or codon == "CAC"):  
            baseCount += 1
    #print (acidCount)
    #print (baseCount)
    if (acidCount == 0 and baseCount == 0):
        print ("There are no acidic or basic codons.")
    elif (acidCount == 0):
        print ("There are no acidic codons, only basic.")
    elif (baseCount == 0):
        print ("There are no basic codons, only acidic.")
    elif (baseCount > acidCount):
        ratio = baseCount / acidCount
        print ("There are more basic codons with a ratio of " + str(ratio) + " basic to 1 acidic.")
    elif (baseCount < acidCount):
        ratio = acidCount / baseCount
        print ("There are more acidic codons with a ratio of " + str(ratio) + " acidic to 1 basic.")
    elif (baseCount == acidCount):
        print ("There are an equal number of acidic and basic codons.")


# Counts the number of times that a sequence of three nucleotides is repeated back to back (eg ATGATG)
def findTriRepeats(DNASeq):
    count = 0
    for i in range (0, len(DNASeq) - 5, 1):
        firstThree = DNASeq[i:i + 3]
        secondThree = DNASeq[i + 3:i + 6]
        if (firstThree == secondThree):
            count += 1
    print ("There are " + str(count) + " times that a three nucleotide sequence is repeated")





if __name__ == "__main__":
    main()
    #end script
    print ("\nend myGeneParser.py")
    exit
