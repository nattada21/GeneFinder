import argparse

parser = argparse.ArgumentParser(description='GeneFinder')
parser.add_argument('sequence', type=str, help='File containing the sequence')
args = parser.parse_args()

def complementary(sequence):
    swap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    transTable = sequence.maketrans(swap)
    complementary_strand = sequence.translate(transTable)
    return complementary_strand

def findORF(seq):
    Long_ORF = ''
    reverse_seq = seq[::-1]
    comp_seq = complementary(seq)
    reverse_comp_seq = comp_seq[::-1]
    allStrands = [seq, reverse_seq, comp_seq, reverse_comp_seq]
    for strand in allStrands:
        ORF = strand[strand.find('ATG'):]
        for i in range(len(ORF)):
            if i%3 == 0: codon = ORF[i:i+3]
            if codon == 'TGA' or codon == 'TAA' or codon == 'TAG':
                ORF = ORF[:i+3]
        if len(ORF) > len(Long_ORF):
            Long_ORF = ORF
    return Long_ORF

def transcribe(ORF):
    mRNA = ORF.replace('T','U')
    return mRNA

def translate(RNA_seq):
    amino_acid = ''
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
    for i in range(len(RNA_seq)):
        if i % 3 == 0:
            codon = RNA_seq[i:i+3]
            for key,value in aa.items():
                if codon in value:
                    amino_acid += key
    return amino_acid

f = open(args.sequence, 'r')
seq = f.read()
#comp = complementary(seq)

ORF = findORF(seq)
mRNA = transcribe(ORF)
print(translate(mRNA))
#print(len('ATGAAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA'))
