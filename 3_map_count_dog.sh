#!/bin/bash

######### FunGen Course Instructions ############
## Purpose: Index the dog reference genome, map trimmed reads with HiSat2, convert and sort alignments with Samtools,
##          count reads with StringTie, and create count matrices with prepDE.py.
## For running on the Alabama Super Computer (ASC): https://hpcdocs.asc.edu/content/slurm-queue-system
## After editing, make executable with "chmod +x map_count_dog.sh" and run with "run_script map_count_dog.sh".
###############################################

########## SLURM Directives
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=map_count_%j.o

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread

########## Define Variables
MyID=aubclsd0338              ## Replace with your ASC ID, e.g., aubtss
WD=/scratch/${MyID}/DogRNAseq   ## Main project folder
CD=${WD}/CleanData              ## Trimmed FASTQ files folder
REFD=${WD}/DogRefGenome         ## Reference genome folder
MAPD=${WD}/Map_HiSat2           ## Mapping output folder
COUNTSD=${WD}/Counts_StringTie  ## StringTie counts folder
RESULTSD=/home/${MyID}/DogRNAseq/Counts_H_S  ## Results in home directory
REF=GCF_000002285.5             ## Tasha the Boxer genome accession

########## Create Directories
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

########## Prepare Reference Genome
###It might have a problem in this step due to the file name( remember)
cd $REFD
cp /home/${MyID}/class_shared/Dog_Tasha_GCF_000002285.5/data/${REF}.fna .
cp /home/${MyID}/class_shared/Dog_Tasha_GCF_000002285.5/data/${REF}.gff .


########## Convert Annotation and Extract Splice Sites/Exons
gffread ${REF}.gff -T -o ${REF}.gtf
hisat2_extract_splice_sites.py ${REF}.gtf > ${REF}.ss
hisat2_extract_exons.py ${REF}.gtf > ${REF}.exon

########## Build HiSat2 Index
hisat2-build --ss ${REF}.ss --exon ${REF}.exon ${REF}.fna Dog_Tasha_index

########## Map Trimmed Reads, Convert, Sort, and Count
cd ${CD}
ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > list
cd ${MAPD}
mv ${CD}/list .
while read i;
do
    ## Map with HiSat2
    hisat2 -p 6 --dta --phred33 \
    -x "${REFD}"/Dog_Tasha_index \
    -1 "${CD}"/"$i"_1_paired.fastq -2 "${CD}"/"$i"_2_paired.fastq \
    -S "$i".sam
    
    ## Convert SAM to BAM and Sort with Samtools
    samtools view -@ 6 -bS "$i".sam > "$i".bam
    samtools sort -@ 6 "$i".bam -o "$i"_sorted.bam
    samtools flagstat "$i"_sorted.bam > "$i"_Stats.txt
    
    ## Count with StringTie
    mkdir "${COUNTSD}"/"$i"
    stringtie -p 6 -e -B -G "${REFD}"/"${REF}".gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i" "${MAPD}"/"$i"_sorted.bam
done < list

########## Create Count Matrices
cd ${COUNTSD}
cp /home/${MyID}/class_shared/prepDE.py3 .
prepDE.py3 -i ${COUNTSD}
cp *.csv ${RESULTSD}
