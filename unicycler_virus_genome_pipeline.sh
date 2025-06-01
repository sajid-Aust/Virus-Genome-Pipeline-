#!/bin/bash
eval "$(conda shell.bash hook)"
## Quality Control
echo "** Viral genome assembly pipeline start **"
echo "** Quality Control Check **"
echo "Do you want to install fastqc (y/n)??"
read choice1

if [[ "$choice1" =~ [yY] ]]; then  
    sudo apt install fastqc
elif [[ "$choice1" =~ [nN] ]]; then  
    echo "Skipping fastqc installation"
else
    echo "Invalid input. Please enter 'y' or 'n'."  # Added else clause
fi

mkdir fastqc_result
echo "Carefully provide the accession number of your fastq file within the 'sra' directory [e.g. SRR1234567]"
read fastq
fastqc -o fastqc_result -t 2 ./sra/"$fastq"*.fastq

echo "* Trimming the fastq files: ***"

echo "Do you want to install trimmomatic (y/n) ??"
read trimmomatic
if [[ "$trimmomatic" =~ [yY] ]]; then  
    sudo apt install trimmomatic
elif [[ "$trimmomatic" =~ [nN] ]]; then  
    echo "Skipping trimmomatic installation"
else
    echo "Invalid input. Please enter 'y' or 'n'."  # Added else clause
fi

mkdir trimmomatic

TrimmomaticPE -phred33 ./sra/"$fastq"_1.fastq ./sra/"$fastq"_2.fastq -baseout trimmomatic/trimmed_"$fastq".fastq LEADING:15 TRAILING:15 SLIDINGWINDOW:4:25 MINLEN:36
echo "Trimming completes succesfully"

echo "# Hybrid Assembly"

# Getting the reference file name
echo "Enter the name of the reference genome file within the 'reference' directory [fasta file] (include the .fasta extension)"
echo "If the reference is an fna file, please rename it to a fasta file (e.g., mv refseq.fna refseq.fasta)"
read refseq

mkdir -p mapped_files # Create the directory if it doesn't exist

# Indexing the Reference Genome
echo "# Indexing the Reference Genome"
bwa index ./reference/"$refseq"

# Aligning Reads to the Reference Genome
echo "# Aligning the trimmed files to the reference genome"
bwa mem -t 190 -M ./reference/"$refseq" ./trimmomatic/trimmed_"$fastq"_1P.fastq ./trimmomatic/trimmed_"$fastq"_2P.fastq > ./mapped_files/bwa.sam


# Extracting Mapped Reads
echo "** Extracting mapped reads **"
samtools view -b -F 4 -o ./mapped_files/mapped.bam ./mapped_files/bwa.sam

# Converting from BAM to FASTQ
echo "** Converting from BAM to FASTQ **"
samtools bam2fq ./mapped_files/mapped.bam > ./mapped_files/mapped.fastq

# Splitting FASTQ file
echo "** Splitting FASTQ file **"
awk 'NR%2==1 { print $0 "/1" } ; NR%2==0 { print substr($0,0,length($0)/2) }' ./mapped_files/mapped.fastq > ./mapped_files/mapped_1.fastq
awk 'NR%2==1 { print $0 "/2" } ; NR%2==0 { print substr($0,length($0)/2+1) }' ./mapped_files/mapped.fastq > ./mapped_files/mapped_2.fastq

mkdir assembly_output
# De Novo Assembly using unicycler
echo "** assembly using unicycler **"
# Activate the 'unicycler' conda environment
source activate unicycler
unicycler -1 ./mapped_files/mapped_1.fastq -2 ./mapped_files/mapped_2.fastq -o assembly_output -t 190

cp ./assembly_output/assembly.fasta ./
# Deactivate the conda environment
conda deactivate 

## Genome annotation
echo "** Genome Annotation **"

# Activate the 'prokka' conda environment
source activate prokka

# Perform genome annotation using Prokka
echo "** Running Prokka annotation (please wait...) **"
prokka --kingdom Viruses --cpus 190 --prefix annotated_genome --locustag annotated_genome scaffolds.fasta

rm ./scaffolds.fasta
# Deactivate the conda environment
conda deactivate  

echo "** Analysis Complete **"

echo "Thank you!"





