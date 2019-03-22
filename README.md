#defines the transcribe function, which turns DNA into mRNA
import matplotlib.pyplot as plt

def Transcribe_DNA(sequence):
    rna_seq = sequence.replace('T','U')
    return rna_seq

def Translate(sequence):
    key = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N",
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S",
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I",

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H",
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R",
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D",
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G",
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V",

                "UAA":"_", "UAC":"Y", "UAG":"_", "UAU":"T",
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S",
                "UGA":"_", "UGC":"C", "UGG":"W", "UGU":"C",
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}

    #Finds the start codon
    start_codon = sequence.find('AUG')

    #Creates a string of the sequence, starting with start codon
    seq_start = sequence[int(start_codon):]

    protein_seq = ''
    sequence_indices = range(len(sequence))
    for n in sequence_indices[0::3]:
        current_codon = sequence[n:n+3]
        if current_codon in key:
            protein_seq += key[current_codon]
    return protein_seq

def Hydrophobicity(protein_seq):
    key1 = {"I":-0.31, "L":-0.56, "V":0.070, "F":-1.13, "M":-0.23, "C":-0.24, "A":0.170,
    "G":0.010, "P":0.450, "T":0.140, "S":0.130, "Y":-0.940, "W":-1.850,
    "Q":0.580, "N":0.420, "H":0.170, "E":2.020, "D":1.230, "K":0.990, "R":0.810}

    seq_start = protein_seq[0:]

    hydro_matrix = ""
    protein_indices = range(len(protein_seq))

    for n in protein_indices[0:]:
        current_aa = str(protein_seq[n])

        if current_aa in key1:

            hydro_matrix += str(key1[current_aa])
            hydro_matrix += " "
    return hydro_matrix

def Matrix(hydro_val):

    hydro_vals = hydro_matrix.split()

    return hydro_vals

def Average(hydro_vals):
    for n in hydro_vals[-1]:
        hydro_floats = map(float, hydro_vals)
        average = sum(hydro_floats)/len(hydro_vals)

    return round(average,2)

#Transcribes and Translates Parsed sequences
#adds them to new variable seq
#prints the protein sequence (seq)

seq = ""

fasta = open('sample2.txt', "r")

seq_ids = []
seqs = []

n = -1

print("\n","Parsing...","\n")
print("\n","Translating...","\n")
print("\n","Transcribing...","\n")

for line in fasta.readlines():

    if line[0] == '>':
        curr_seq_id = line[1:]
        seq_ids.append(curr_seq_id)
        seqs.append('')
        n += 1
    else:
        seqs[n] += line
        rna = Transcribe_DNA(line)
        protein = Translate(rna)
        seq += (protein)

        continue
#calls all functions for each sequence
#displays protein sequence and hydrophobicity data
for seq_id, seq in zip(seq_ids, seqs):
    print("**Protein Sequence for:", seq_id)
    rna = Transcribe_DNA(seq)
    protein = Translate(rna)
    hydro_matrix = Hydrophobicity(protein)
    hydro_vals = Matrix(hydro_matrix)
    average = Average(hydro_vals)
    plt.plot(average)
    plt.ylabel("Hydrophobicity")
    plt.show()
    print("\t", protein, "\n")
    print("**Hydrophobicity data for:", seq_id)
    print("\t","Average interface scale value:", average, "\n")
    print("Matrix:", "\n","\t", hydro_vals, "\n")
