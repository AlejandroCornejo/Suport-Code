

"""
Este fichero se emplea para leer el archivo de intersecciones y escribe los bloques de elementos de acuerdo con las familias H, V y DF
"""

input_name = "Asco_I.intersections"
output_name = "bloques_elementos.txt"

input_file  = open(input_name, "r")
output_file = open(output_name, "w")

lines = input_file.readlines()

horizontal_elems = []
vertical_elems = []
cupula_elems = []

print("Leyendo fichero de intersecciones...")

intersection_block = False
for line in lines:
    split_line = line.split()
    if len(split_line) > 0: # Skip empty lines
        if line.find("intersection:") != -1: # a new tendon 
            intersection_block = True
            tendon_name = split_line[5]
        elif line.find("End") != -1 and line.find("intersection") != -1:
            intersection_block = False
        if intersection_block and line.find("Begin") == -1: # Id's bloque
            id_elem = int(split_line[0])
            if tendon_name[7] == "H":
                horizontal_elems.append(id_elem)
            elif tendon_name[7] == "V":
                vertical_elems.append(id_elem)
            elif tendon_name[7] == "D":
                cupula_elems.append(id_elem)
print("Escribiendo elementos...")
output_file.write("\n \n Elementos Tendones Horizontales: \n \n")

colum = 0
for elem in horizontal_elems:
    output_file.write("  " + str(elem))
    colum += 1
    if colum == 30:
        output_file.write("\n")
        colum = 0

output_file.write("\n \n Elementos Tendones Verticales: \n \n")
colum = 0
for elem in vertical_elems:
    output_file.write("  " + str(elem))
    colum += 1
    if colum == 30:
        output_file.write("\n")
        colum = 0

output_file.write("\n \n Elementos Tendones Cupula: \n \n")

colum = 0
for elem in cupula_elems:
    output_file.write("  " + str(elem))
    colum += 1
    if colum == 30:
        output_file.write("\n")
        colum = 0

input_file.close()
output_file.close()


