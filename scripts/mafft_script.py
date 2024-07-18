import os
import subprocess

for file in os.listdir('/multifasta'):
    gene = file.replace('multi_', '.fasta', '')
    command = f"'/cluster/software/MAFFT/7.520-GCC-12.3.0-with-extensions/bin/mafft' --auto --clustalout --reorder '/multifasta/{file}' > 'alignments/align_{gene}.MSA'"
    command = f"/cluster/software/MAFFT/7.520-GCC-12.3.0-with-extensions/bin/mafft  --auto --clustalout --reorder "{gene}" > "/alignments/align_{gene}.MSA""
    print(f'Alignment of {gene} correctly done')
    result = subprocess.run(command, shell=True, capture_output=True, text=True)