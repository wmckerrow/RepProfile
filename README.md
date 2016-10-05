# RepProfile
RepProfile takes advantages of the depenedencies between read seqeunces to provide a more accurate alignment of paired end reads to repetitive sequence. If you use RepProfile in your research please cite our paper: <PAPER INFO>.

For an overview of the method used by RepProfile, see our paper:
<PAPER INFO>
For a more detailed description of the mathematics used by RepProfile, see the paper supplement.

# Installation
RepProfile requires no compiling, but has the following python library dependencies:
BIO
pysam
numpy
scipy

The QuickStart example also requires the following be installed:
bwa
samtools

# Quick Start
To get started quickly, align the QuickStart example as follows. Edit the first three lines of make_candidate_alignments.sh to use your own read data and reference. RepProfile is not fast enough for genomewide alignment, so use a reference that only includes one repeat type along with flanking sequence that excedes that length of almost all fragments. See manual.md for details.

```bash
utilities/make_candidate_alignments.sh 
python main/parse_PEbam_for_RepProfile.py -b candidates.masked.aln.bam -r QuickStart/reads_R1.fastq,QuickStart/reads_R2.fastq -g QuickStart/someFB.fa -m 4
python main/RepProfile.py -a alignments_list.txt -r QuickStart/someFB.fa -p QuickStart/HyperEditingPrior.txt -n 5
```

You can generate a bam file for the alignment:
samtools view -H candidates.masked.aln.bam > QS.sam
python utilities/write_sam.py >> QS.sam
samtools view -bS QS.sam | samtools sort - -o QS.bam
rm QS.sam

Or report sequences/positions predicted to be hyper-edited:
python utilities/report.py -r hyper_f,hyper_r
python utilities/report.py -p edit_f,edit_r

You can clean up with:
rm alignments_*.pkl alignments_list.txt candidates.masked.aln.bam expU_f.pkl expU_r.pkl f_prob.pkl pos_type.pkl rep_type.pkl r_prob.pkl genome_profile_f.pkl genome_profile_initial.pkl genome_profile_r.pkl

More detailed information regarding how to run RepProfile can be found in manual.md.
