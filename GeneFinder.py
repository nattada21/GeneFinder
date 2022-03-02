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
    reverse_seq = seq[::-1]
    comp_seq = complementary(seq)
    reverse_comp_seq = comp_seq[::-1]
    
    ORF = ''
    '''if 'ATG' in strand:
        ORF = seq[seq.find('ATG'):]
        Stop_1 = ORF_seq.find('TAA')
        Stop_2 = ORF_seq.find('TAG')
        Stop_3 = ORF_seq.find('TGA')'''

    return ORF

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
#print(seq)
#print(comp)

ORF = findORF(seq)
mRNA = transcribe(ORF)
print(translate(mRNA))
