#!/usr/bin/env python

#from intermine.webservice import Service
import subprocess
import os
import shutil
import re
module = 'module load BLAST+/2.14.1-gompi-2023a'
subprocess.run(module, shell=True, capture_output=True, text=True)

"""
print('\nSTABLISHING CONECTION WITH YEASTMINE...')
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")
genes = list()

with open('genes_2.txt', 'r') as f:
    print("\nREADING FILE WITH THE GENES' NAMES...")
    for line in f:
        genes.append(line.strip('\n'))

other_genes = list()
for gene in os.listdir():
    if gene.endswith('.fasta') and not gene.startswith('maltose'):
        other_genes.append(gene)        

with open('maltose_genes.fa', 'w') as f:
    print('\nSEARCHING SEQUENCES...')
    for gene in genes:
        print(f'\nSEARCHING {gene}')
        query = service.new_query("Gene")
        query.add_view(
            "primaryIdentifier", "secondaryIdentifier", "symbol", "name",
            "organism.shortName", "sequence.length", "sequence.residues"
        )
        query.add_constraint("Gene", "LOOKUP", gene, code="E")
        query.add_constraint("status", "=", "Active", code="C")
        query.add_constraint("status", "IS NULL", code="D")
        query.set_logic("(C or D) and E")
        for row in query.rows():
            print(f'\nWRITING {gene} GENE IN FASTA FILE...')
            f.write(f">{row['symbol']}\n{row['sequence.residues']}\n")
    
    for gene in other_genes:
        print(f'\nWRITING {gene} GENE IN FASTA FILE...')
        with open(gene, 'r') as g:
            for line in g:
                if line.startswith('>'):
                    f.write('\n' + line)
                else:
                    f.write(line.strip('\n'))

files = os.listdir('genomes/')
seqs = ''
print('\n')
for file in files:
    if file.endswith('.fasta'):
        with open(f'genomes/{file}', 'r') as a:
            print(f'Reading {file}...')
            for line in a:
                seqs += line
    print('\n')
"""
"""with open('compilation.fasta', 'w') as c:
    print('\nCREATING MULTIFASTA...')
    c.write(seqs)"""


# Diccionario para convertir números a números romanos
def int_to_roman(num):
    val = [
        1000, 900, 500, 400,
        100, 90, 50, 40,
        10, 9, 5, 4,
        1
        ]
    syb = [
        "M", "CM", "D", "CD",
        "C", "XC", "L", "XL",
        "X", "IX", "V", "IV",
        "I"
        ]
    roman_num = ''
    i = 0
    while  num > 0:
        for _ in range(num // val[i]):
            roman_num += syb[i]
            num -= val[i]
        i += 1
    return roman_num

# Función para reemplazar 'chr' seguido de un número con 'chr' seguido del número romano
def replace_chr_with_roman(content):
    def replacer(match):
        num = int(match.group(1))
        return 'chr' + int_to_roman(num)

    return re.sub(r'chr(\d+)', replacer, content)

# Leer el fichero, reemplazar y escribir de nuevo
def replace_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    new_content = replace_chr_with_roman(content)

    with open(file_path, 'w') as file:
        file.write(new_content)




"""print('\nCREATING DATABASE...')
db = "makeblastdb -in compilation.fasta -dbtype nucl metadata_output_prefix db.njs -out compilation_db"
subprocess.run(db, shell=True, capture_output=True, text=True)
"""
with open("header.txt", "w") as f:
    header = "gene\taccession\t%match\talign_length\tmismatch\tgapopen\tquery_length\tstart\tend\te-value\tbitscore\n"
    f.write(header)


with open('genes.fa', 'r') as f:
    n = list()
    with open('genes.gff3', 'w') as g:
        gff_header =  'seqid\tsource\ttype\tstart\tend\tlength\tstrand\tphase\tattributes\n'
        g.write(gff_header)

    for line in f:
        if line.startswith('>'):
            name = line[1:].strip('\n')
            name = name.lstrip('>')
            name = name.split(' ')
            name = name[0]
            n.append(name)
            gene = ''
            gene += line
        else:
            gene += line.strip('\n')

            with open(f'{name}.fa', 'w') as g:
                g.write(gene)



for name in n:
    print(f'\nDOING BLAST OF {name}...')
    command = f'blastn -db compilation_db -query {name}.fa -out  ~/TFM/blast/results_{name}.out -perc_identity 95 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    command_gff = f'python3 blast2gff.py -b  ~/TFM/blast/results_{name}.out >>  genes.gff3'
    result = subprocess.run(command_gff, shell=True, capture_output=True, text=True)
    """with open (f'results_{name}.out', 'r') as o:
                lines = o.readlines()
                for line in lines[1:]:
                    cols = line.strip().split('\t')
                    accession = cols[1]
                    start = cols[8]
                    end = cols[9]
                    if start <= end:
                        strand = 'plus'
                    else:
                        strand = 'minus'
                    #print(accession)
                    command_fasta = f"blastdbcmd -db compliation_db -entry {accession} -range {start}-{end} -strand {strand} -outfmt %f -out {name}_hits.fa"
                    print(command_fasta)
                    subprocess.run(command_fasta, shell=True)
                    with open(f'{name}_hits.fa', 'r') as h:
                        with open('multifasta_genes.fa', 'a') as multifasta:
                            shutil.copyfileobj(h, multifasta)"""

    if result.returncode == 0:
        print(f"BLASTN for {name}.fa was correctly executed.")
    else:
        print(f"Error in BLASTN for {name}.fa:")
        print(result.stderr)

blast = list()
for i in n:
    subprocess.run(f"cat header.txt  ~/TFM/blast/results_{i}.out > ~/TFM/blast/{i}_blast.tsv", shell=True)
    subprocess.run(f"rm  ~/TFM/blast/results_{i}.out", shell=True)
    blast.append(f'/cluster/home/lauvapo/TFM/blast/{i}_blast.tsv')

subprocess.run("rm header.txt", shell=True)

# Usar la función con el nombre del fichero
for file in blast:
    print('Sustituyendo...')
    replace_in_file(file)
    print(f'{file} modificado')
