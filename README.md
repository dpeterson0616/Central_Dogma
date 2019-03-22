
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

