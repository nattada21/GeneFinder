import argparse

#Uses argparse to accept user input
#User needs to enter a file name containing the sequence
parser = argparse.ArgumentParser(description='GeneFinder')
parser.add_argument('sequence', type=str, help='File containing the sequence')
args = parser.parse_args()

#Finds the complementary DNA strand given a DNA strand
def complementary(sequence):
    #Dictionary indicating the nucleotides and their corresponding base pairs
    swap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    transTable = sequence.maketrans(swap)
    #Creates a complementary strand with the swaps
    complementary_strand = sequence.translate(transTable)
    return complementary_strand

#Finds the longest open reading frame
def findORF(seq):
    Long_ORF = '' #Empty string that will hold longest ORF later
    reverse_seq = seq[::-1] #Reverse string from the original DNA strand
    comp_seq = complementary(seq) #Calls complementary function for the complementary string
    reverse_comp_seq = comp_seq[::-1] #Reverse of the complementary strand
    allStrands = [seq, reverse_seq, comp_seq, reverse_comp_seq] #creates a list with all the strands
    for strand in allStrands: #Iterates through each strand in the list
        ORF = strand[strand.find('ATG'):] #Splits the DNA strand where the start codon is located
        for i in range(len(ORF)): #Iterates through ORF to create codons
            if i%3 == 0: codon = ORF[i:i+3] #Creates codons
            if codon == 'TGA' or codon == 'TAA' or codon == 'TAG': #Checks if codon is a stop codon
                ORF = ORF[:i+3] #Assigns ORF to the substring from start to stop codons
        if len(ORF) > len(Long_ORF): #compares the ORF to the Longest ORF
            Long_ORF = ORF #Assigns the longest ORF variable
    return Long_ORF

#Transcribes the longest ORF into mRNA sequence
def transcribe(ORF):
    mRNA = ORF.replace('T','U')
    return mRNA

#Translates mRNA sequence into a protein sequence
def translate(RNA_seq):
    amino_acid = '' #Empty string that will hold the amino acids
    #Library containing all the amino acids and their corresponding codons
    aa = {'A':['GCU','GCC','GCA','GCG'],
          'R':['CGU','CGC','CGA','CGG','AGA','AGG'],
          'N':['AAU','AAC'],
          'D':['GAU','GAC'],
          'C':['UGU','UGC'],
          'E':['GAA','GAG'],
          'Q':['CAA','CAG'],
          'G':['GGU','GGC','GGA','GGG'],
          'H':['CAU','CAC'],
          'I':['AUU','AUC','AUA'],
          'L':['UUA','UUC','CUU','CUC','CUA','CUG'],
          'K':['AAA','AAG'],
          'M':['AUG'],
          'F':['UUU','UUC'],
          'P':['CCU','CCC','CCA','CCG'],
          'S':['UCU','UCC','UCA','UCG','AGU','AGC'],
          'T':['ACU','ACC','ACA','ACG'],
          'W':['UGG'],
          'Y':['UAU','UAC'],
          'V':['GUU','GUC','GUA','GUG']}
    for i in range(len(RNA_seq)): #Goes through the length of the mRNA sequence
        if i % 3 == 0:
            codon = RNA_seq[i:i+3] # Creates codons
            for key,value in aa.items(): #Separates keys and values
                if codon in value: #Checks to see if a codon is in the dictionary
                    amino_acid += key #Adds the corresponding amino acid to the string
    return amino_acid

#Opens file containing DNA strand from 5' to 3'
f = open(args.sequence, 'r')
#Reads file and assigns contents to a variable
seq = f.read()
#Finds the longest Open reading frame
ORF = findORF(seq)
#Produces and mRNA sequnce given the longest open reading frame
mRNA = transcribe(ORF)
#Prints the protein sequence given the mRNA strand produced above
print(translate(mRNA))
