export reads1=QuickStart/reads_R1.fastq #Reads that (might) align to the repeat of interest. See manual.md for details.
export reads2=QuickStart/reads_R2.fastq #Reads that (might) align to the repeat of interest.
export ref=QuickStart/someFB.fa #Repeat genome. One fasta entry per repeat instance. Each entry should include the repeat sequence with (1kb) of flanking seqeunce upstream and downstream.

python utilities/mask_fq.py $reads1 T C > $reads1.TCmasked
python utilities/mask_fq.py $reads2 A G > $reads2.AGmasked
python utilities/mask_fa.py $ref A G > $ref.AGmasked
python utilities/mask_fa.py $ref T C > $ref.TCmasked
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
