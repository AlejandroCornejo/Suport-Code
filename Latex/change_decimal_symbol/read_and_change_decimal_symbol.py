'''
Este script lee un fichero y cambia el simbolo decimal
a todos aquellos números 
'''

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


nombre_input  = "input.tex"
nombre_output = "output.tex"

input_file  = open(nombre_input, "r")
output_file = open(nombre_output, "w")

old_decimal_symbol = '.'
new_decimal_symbol = ','

# leemos lineas
lines = input_file.readlines()
line_counter = 0
# bucle sobre las lineas del doc
for line in lines:
    line_counter += 1
    # separamos la linea en palabras
    split_line = line.split()
    # saltamos lineas vacias
    if len(split_line) > 0:
        for word in split_line:
            if is_number(word):
                print("Number changed: ", word, " in line ", line_counter)
                word = word.replace(old_decimal_symbol, new_decimal_symbol)
            output_file.write(word + " ")
        output_file.write("\n")
    else:
        output_file.write("\n")

