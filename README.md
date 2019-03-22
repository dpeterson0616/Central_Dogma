#defines the transcribe function, which turns DNA into mRNA

def Transcribe_DNA(sequence):
    rna_seq = sequence.replace('T','U')
    print(rna_seq)
    return rna_seq

