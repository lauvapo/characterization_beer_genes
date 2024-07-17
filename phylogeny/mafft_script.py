import os
import subprocess

for file in os.listdir('multifasta'):
    gene = file.replace('.fasta', '')
    command = (f"mafft --thread 12 --globalpair --maxiterate 1000 --clustalout --reorder 'multifasta/{file}' > 'align_{gene}.MSA'")
    print(f'Alignment of {gene} correctly done')
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
