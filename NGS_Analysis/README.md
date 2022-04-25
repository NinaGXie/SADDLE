Note that this NGS Analysis folder contains a copy of the Matlab and python code of the SADDLE NGS data analysis pipeline. 
The MATLAB code used for multiplex PCR primer algorithm is available upon request to the corresponding author under NDA for an academic lab.

1. Run Data_process.m to do adapter trimming and quality control of the raw NGS .fastq file. In this step, we sort all the reads into long reads (>60) and short reads(<60). Long reads file will be used for Step2 alignment.
2. Run bowtie2 to do alignment, we don’t provide bowtie2 here. Please 1. Build your own reference, and 2. Run alignment to generate a .sam file. An alignment example can be ‘bowtie2 -x Ref -1 file_F.fastq -2 file_R.fasta -S output.sam’
3. Run analysis_code.ipynb to generate dimer candidates file for further analysis in step 4
4. Run Dimer_Identify.m to generate dimer matrix.