export reads1=/users/wmckerro/data/wmckerro/YAS_reads/PE_RNA_genewiz/4_GAATTCGT-GGCTCTGA_L002_R1_001.fastq.gz #Illumina reads. Gunzipped.
export reads2=/users/wmckerro/data/wmckerro/YAS_reads/PE_RNA_genewiz/4_GAATTCGT-GGCTCTGA_L002_R2_001.fastq.gz #Illumina reads. Gunzipped.
export ref=/users/wmckerro/data/wmckerro/dm6/dm6.fa #Reference genome.
export repeatbed=/users/wmckerro/data/wmckerro/fly_annotations/dm6.reps.bed #List of repeat instances in bed format. Can download from UCSC genome table browser.
export repeatname=FB4_DM #Name of repeat to generate alignments for.
export RepProfile=/users/wmckerro/RepProfile

#module load bwa
#module load samtools

python $RepProfile/utilities/make_repeat_genome.py -r $repeatname -f 1000 -m $repeatbed -g $ref > repeatgenome.fa

python $RepProfile/utilities/mask_fq_fromgz.py $reads1 T C > $reads1.TCmasked
python $RepProfile/utilities/mask_fq_fromgz.py $reads2 A G > $reads2.AGmasked
python $RepProfile/utilities/mask_fa.py $ref A G > $ref.AGmasked
python $RepProfile/utilities/mask_fa.py $ref T C > $ref.TCmasked
bwa index $ref.AGmasked
bwa index $ref.TCmasked
bwa aln -t 8 -n 4 -l 25 -i 10 $ref.AGmasked $reads1.TCmasked >  AG.masked.R1.aln.sai
bwa aln -t 8 -n 4 -l 25 -i 10 $ref.AGmasked $reads2.AGmasked > AG.masked.R2.aln.sai
bwa sampe $ref.AGmasked AG.masked.R1.aln.sai  AG.masked.R2.aln.sai $reads1.TCmasked $reads2.AGmasked > AG.masked.aln.sam
samtools view -bS AG.masked.aln.sam > AG.masked.aln.bam
bwa aln -t 8 -n 4 -l 25 -i 10 $ref.TCmasked $reads1.TCmasked > TC.masked.R1.aln.sai
bwa aln -t 8 -n 4 -l 25 -i 10 $ref.TCmasked $reads2.AGmasked > TC.masked.R2.aln.sai
bwa sampe $ref.TCmasked TC.masked.R1.aln.sai TC.masked.R2.aln.sai $reads1.TCmasked $reads2.AGmasked > TC.masked.aln.sam
samtools view -bS TC.masked.aln.sam > TC.masked.aln.bam
samtools merge -f masked.aln.unsort.bam AG.masked.aln.bam TC.masked.aln.bam
samtools sort masked.aln.unsort.bam > masked.aln.bam
samtools index masked.aln.bam

rm $reads1.TCmasked $reads2.AGmasked
rm AG.masked.R1.aln.sai AG.masked.R2.aln.sai AG.masked.aln.sam AG.masked.aln.bam
rm TC.masked.R1.aln.sai TC.masked.R2.aln.sai TC.masked.aln.sam TC.masked.aln.bam
rm masked.aln.unsort.bam

python $RepProfile/utilities/repeat_read_ids.py -r $repeatbed -b masked.aln.bam -i $repeatname > repeat_read_ids.txt
python $RepProfile/utilities/make_repeat_fastq.py -l repeat_read_ids.txt -g $reads1 > reads_R1.fastq
python $RepProfile/utilities/make_repeat_fastq.py -l repeat_read_ids.txt -g $reads2 > reads_R2.fastq

export reads1=reads_R1.fastq
export reads2=reads_R2.fastq
export ref=repeatgenome.fa

python $RepProfile/utilities/mask_fq.py $reads1 T C > $reads1.TCmasked
python $RepProfile/utilities/mask_fq.py $reads2 A G > $reads2.AGmasked
python $RepProfile/utilities/mask_fa.py $ref A G > $ref.AGmasked
python $RepProfile/utilities/mask_fa.py $ref T C > $ref.TCmasked
bwa index $ref.AGmasked
bwa index $ref.TCmasked
bwa aln -t 8 -N -n 4 -l 25 -i 10 $ref.AGmasked $reads1.TCmasked >  candidates.AG.masked.R1.aln.sai
bwa aln -t 8 -N -n 4 -l 25 -i 10 $ref.AGmasked $reads2.AGmasked > candidates.AG.masked.R2.aln.sai
bwa sampe -n 1000 -N 0 $ref.AGmasked candidates.AG.masked.R1.aln.sai  candidates.AG.masked.R2.aln.sai $reads1.TCmasked $reads2.AGmasked > candidates.AG.masked.aln.sam
samtools view -bS candidates.AG.masked.aln.sam > candidates.AG.masked.aln.bam
bwa aln -t 8 -N -n 4 -l 25 -i 10 $ref.TCmasked $reads1.TCmasked > candidates.TC.masked.R1.aln.sai
bwa aln -t 8 -N -n 4 -l 25 -i 10 $ref.TCmasked $reads2.AGmasked > candidates.TC.masked.R2.aln.sai
bwa sampe -n 1000 -N 0 $ref.TCmasked candidates.TC.masked.R1.aln.sai candidates.TC.masked.R2.aln.sai $reads1.TCmasked $reads2.AGmasked > candidates.TC.masked.aln.sam
samtools view -bS candidates.TC.masked.aln.sam > candidates.TC.masked.aln.bam
samtools merge -f candidates.masked.aln.unsort.bam candidates.AG.masked.aln.bam candidates.TC.masked.aln.bam
samtools sort -n candidates.masked.aln.unsort.bam > candidates.masked.aln.bam

rm $reads1.TCmasked $reads2.AGmasked
rm $ref.AGmasked $ref.AGmasked.amb $ref.AGmasked.ann $ref.AGmasked.bwt $ref.AGmasked.pac $ref.AGmasked.sa
rm $ref.TCmasked $ref.TCmasked.amb $ref.TCmasked.ann $ref.TCmasked.bwt $ref.TCmasked.pac $ref.TCmasked.sa
rm candidates.AG.masked.R1.aln.sai candidates.AG.masked.R2.aln.sai candidates.AG.masked.aln.sam
rm candidates.TC.masked.R1.aln.sai candidates.TC.masked.R2.aln.sai candidates.TC.masked.aln.sam
rm candidates.AG.masked.aln.bam candidates.TC.masked.aln.bam candidates.masked.aln.unsort.bam
