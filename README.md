def Matrix(hydro_val):

    hydro_vals = hydro_matrix.split()

    return hydro_vals

def Average(hydro_vals):
    for n in hydro_vals[-1]:
        hydro_floats = map(float, hydro_vals)
        average = sum(hydro_floats)/len(hydro_vals)

    return round(average,2)

