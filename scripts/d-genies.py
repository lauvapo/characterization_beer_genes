# PROCESA DATOS DE D-GENIES
# LE PASAMOS LOS ARCHIVOS QUE CONTIENEN LAS ASOCIACIONES ENTRE MULTIFASTA DE CADA GEN CON GENOMA DE REFERENCIA
# SUSTITUIMOS EL NÚMERO DE ACCESO DEL S288C POR EL CROMOSOMA QUE CORRESPONDA CON EL ARCHIVO sequence_report.tsv
# GENERA TSV PARA CADA GEN PONIENDO EL CROMOSOMA DONDE HA DADO MATCH EL CONTIG FRENTE A LA CEPA S288C
import os

# Directorio y archivos
dir = 'C:/Users/Usuario/Downloads/'
dir2= 'C:/Users/Usuario/Downloads/asociaciones'
align = 'genes_GCA_000146045.2_R64_genomic_assoc(1).tsv'
ref = 'sequence_report.tsv'
output_file = 'aligned_output.tsv'
'_multi_GCA_000146045.2_R64_genomic_assoc'
chromosomes = dict()

# Leer el archivo de referencia y construir el diccionario de cromosomas
with open(os.path.join(dir, ref), 'r') as r:
    for line in r:
        line = line.strip().split('\t')
        if line[6].startswith('BK0069'):
            chromosomes[line[6]] = line[12]

# Leer el archivo de alineación y escribir el archivo de salida con las sustituciones
with open(os.path.join(dir, align), 'r') as a, open(os.path.join(dir, output_file), 'w') as out:
    for line in a:
        line_split = line.strip().split('\t')
        key = line_split[1]
        if key in chromosomes:
            line_split[1] = chromosomes[key]
        out.write('\t'.join(line_split) + '\n')

print(f"Archivo procesado y guardado en {os.path.join(dir, output_file)}")


for file in os.listdir(dir2):
    gen = file.replace('_multi_GCA_000146045.2_R64_genomic_assoc(1).tsv', '')
    with open(os.path.join(dir2,file), 'r') as a, open(os.path.join(dir2, f'{gen}_chr.tsv'), 'w') as out:
        for line in a:
            line_split = line.strip().split('\t')
            key = line_split[1]
            if key in chromosomes:
                line_split[1] = chromosomes[key]
            out.write('\t'.join(line_split) + '\n')



