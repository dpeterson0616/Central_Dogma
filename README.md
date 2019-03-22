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
