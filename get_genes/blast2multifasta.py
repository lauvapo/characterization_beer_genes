# Script for converting blast hits into multifasta sequences for each gene

import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Directorios
dir_blast = '/cluster/home/lauvapo/TFM/blast/'
dir_genomes = '/cluster/home/lauvapo/TFM/genomes/'
output_dir = '/cluster/home/lauvapo/TFM/multifasta/'

# Crear el directorio de salida si no existe
os.makedirs(output_dir, exist_ok=True)

# Listar archivos
blast_files = os.listdir(dir_blast)
genome_files = os.listdir(dir_genomes)

# Inicializar diccionario de genes
genes = dict()

# Construir el diccionario de genes a partir de los archivos .tsv
for file in blast_files:
    if file.endswith('.tsv') and not file.startswith('CNV'):
        with open(os.path.join(dir_blast, file), 'r') as f:
            gene = file.replace('_blast.tsv', '')
            for line in f:
                if not line.startswith('accession'):
                    line = line.strip().split('\t')
                    if gene not in genes:
                        genes[gene] = []
                    genes[gene].append([line[0], int(line[1]), int(line[2])])


# Inicializar diccionario de coordenadas de contigs
coord_chr = dict()

# Procesar cada gen y sus contigs
for gene, contigs in genes.items():
    for contig_info in contigs:
        contig_name = contig_info[0]
        start = int(contig_info[1])
        end = int(contig_info[2])
        if start > end: # Se hará reverse complement después
            start2 = start
            end2 = end
            end = start2
            start = end2
            
        print(contig_info[0])
        print(f'{start}- {end}')
        # Extraer la cepa del nombre del contig
        strain = re.sub(r'_chr.*', '', contig_name)

        # Saltar ciertas cepas
        if strain.startswith('LIUU01') or strain.startswith('CABIK'): ##### BUSCAR
            continue

        # Leer el archivo fasta correspondiente a la cepa
        fasta_file = os.path.join(dir_genomes, f"{strain}.fasta")
        if os.path.exists(fasta_file):
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                if f'{strain}_{seq_record.id}' == contig_name:
                    seq = seq_record.seq[start-1:end]
                    print(start, end, len(seq))
                    if gene not in coord_chr:
                        coord_chr[gene] = []
                    if not seq.startswith('ATG') and seq.endswith('CAT'):
                        seq = seq.reverse_complement()
                        
                    if len(seq) > 0:
                        coord_chr[gene].append((contig_name, seq))
                    break

# Función para escribir secuencias en un archivo multifasta
def write_multifasta(gene, sequences):
    filename = f"{gene}_multi.fasta"
    filepath = os.path.join(output_dir, filename)
    seq_records = [SeqRecord(seq, id=contig, description="") for contig, seq in sequences]
    with open(filepath, "w") as f:
        if len(seq) > 1:
            SeqIO.write(seq_records, f, "fasta")

# Recorrer el diccionario y escribir cada grupo de secuencias de un gen en un archivo multifasta
for gene, sequences in coord_chr.items():
    write_multifasta(gene, sequences)

print("Archivos multifasta creados en", output_dir)
