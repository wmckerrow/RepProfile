# RepProfile
RepProfile takes advantages of the depenedencies between read seqeunces to provide a more accurate alignment of paired end reads to repetitive sequence. If you use RepProfile in your research please cite our paper: <PAPER INFO>.

For an overview of the method used by RepProfile, see our paper:
<PAPER INFO>
For a more detailed description of the mathematics used by RepProfile, see the paper supplement.

## Installation
RepProfile requires no compiling, but has the following python library dependencies:
BIO
pysam
numpy
scipy

The QuickStart example also requires the following be installed:
bwa
samtools

## QuickStart Example
To make sure that RepProfile is working, align the QuickStart example as follows. Edit the first three lines of make_candidate_alignments.sh to use your own read data and reference. RepProfile is not fast enough for genomewide alignment, so use a reference that only includes one repeat type along with flanking sequence that excedes that length of almost all fragments. See manual.md for details about how to generate this.

```bash
utilities/make_candidate_alignments.sh 
python main/parse_PEbam_for_RepProfile.py -b candidates.masked.aln.bam -r QuickStart/reads_R1.fastq,QuickStart/reads_R2.fastq -g QuickStart/someFB.fa -m 4
python main/RepProfile.py -a alignments_list.txt -r QuickStart/someFB.fa -p QuickStart/HyperEditingPrior.txt -n 5
```

You can generate a bam file for the alignment:
```bash
samtools view -H candidates.masked.aln.bam > QS.sam
python utilities/write_sam.py >> QS.sam
samtools view -bS QS.sam | samtools sort - -o QS.bam
rm QS.sam
```

Or report sequences/positions predicted to be hyper-edited:
```bash
python utilities/report.py -r hyper_f,hyper_r
python utilities/report.py -p edit_f,edit_r
```

You can clean up with:
```bash
rm alignments_*.pkl alignments_list.txt candidates.masked.aln.bam expU_f.pkl expU_r.pkl f_prob.pkl pos_type.pkl rep_type.pkl r_prob.pkl genome_profile_f.pkl genome_profile_initial.pkl genome_profile_r.pkl
```

More detailed information regarding how to run RepProfile can be found in manual.md.

## Generate repeat genome and reads
Given a named repeat, a bed file with repeats (can be generated on the UCSC genome table browser,) and a reference genome, you can generate a repeat genome as follows:
```bash
python utilities/make_repeat_genome.py -r FB4_DM -f 1000 -m dm6.reps.bed -g dm6.fa > repeatgenome.fa
```
We recommend that you use only reads that are known to align to the repeat of interest.Given an alignment in bam format, you can collect reads that align to the given repeat as follows:
```bash
python utilities/repeat_read_ids.py -r dm6.reps.bed -b alignment.bam -i FB4_DM > repeat_read_ids.txt
python utilities/make_repeat_fastq.py -l repeat_read_ids.txt -g illumina_reads_R1.fastq.gz > reads_R1.fastq
python utilities/make_repeat_fastq.py -l repeat_read_ids.txt -g illumina_reads_R2.fastq.gz > reads_R2.fastq
```

If you are looking for hyper editing, utilities/make_candidate_alignments_genomic.sh will automatically build a repeat genome, find repeat reads, and build candidate alignments. Edit the first four lines appropriately.

## Contact
Questions? email Will McKerrow: wilson_mckerrow@brown.edu.