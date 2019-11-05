module load PEAR

mkdir data/bylib/fqt_merge

pear -f GARS_cDNA.read1.fq.gz \
     -r GARS_cDNA.read2.fq.gz -v 20 -j 2 \
     -o stephanie_GARS_cDNA

# map assembled reads to junctions
echo "\
>wt_junc
GAAATGATCTATCCCCTCCAGTGTCTTTTAACTTAATGTTCAAGACTTTCATTGGGCCTGGAGGAAACATGCCTGGGTACTTGAGACCAGAAACTGCACAGGGGATTTTCTTGAATTTCAAACGACTTTTGGAGTTCAACCAAGGAAAGTTGCCTTTTGCTGCTGCCCAGATTGGAAATTCTTTTAGAAATGAGATCTCCCCTCGATCTGGACTGATCAGAGT
>del_junc
GAAATGATCTATCCCCTCCAGTGTCTTTTAACTTAATGTTCAAGACTTTCATTGGGCCTGGAGGAAACATGCCTGGGTACTTGAGACCGGGGATTTTCTTGAATTTCAAACGACTTTTGGAGTTCAACCAAGGAAAGTTGCCTTTTGCTGCTGCCCAGATTGGAAATTCTTTTAGAAATGAGATCTCCCCTCGATCTGGACTGATCAGAGT\
">junction_ref.fa

head -n2 junction_ref.fa > wt_junc.fa
bwa index -a is wt_junc.fa
tail -n2 junction_ref.fa > del_junc.fa
bwa index -a is del_junc.fa

# map reads against each individual junction reference

bwa mem -t 30 wt_junc.fa stephanie_GARS_cDNA.assembled.fastq | samtools view -bS -F 2048 - > merged_vs_wtjunc.bam
bwa mem -t 30 del_junc.fa stephanie_GARS_cDNA.assembled.fastq | samtools view -bS -F 2048  - > merged_vs_deljunc.bam

# compare alignment status, map quality, and alignment score to determine which junction to assign to

python compare_us_bams_two_refs.py \
  --inBam1 merged_vs_wtjunc.bam \
  --inBam2 merged_vs_deljunc.bam \
  --outBam1 merged_vs_wtjunc.an.bam \
  --outBam2 merged_vs_deljunc.an.bam \
  --mapqMinThresh 0 --mapqDiffThresh 1 \
  --flagRequireEither 0 \
  --useASinsteadofNM --minNumMismatchDiffThresh 1 \
  --outSingletonOnlySummary out.1 \
  --outMatedOnlySummary  out.2 \
  --minNumMismatchDiffThresh 2 \
  --outSummary out.3 \
  --refname1 wt --refname2 del 

# output from resulting report

"""
822423  isMappedBoth    mqwtpass        mqdelpass       mqWithinThresh  nmWithinThresh
689502  isMappedBoth    mqwtpass        mqdelpass       mqWithinThresh  nmwtdiffoverThresh
595223  isMappedBoth    mqwtpass        mqdelpass       mqWithinThresh  nmdeldiffoverThresh
33068   isMappedBoth    mqwtpass        mqdelpass       mqdeldiffoverThresh     nmdeldiffoverThresh
32436   isMappedBoth    mqwtpass        mqdelpass       mqwtdiffoverThresh      nmwtdiffoverThresh
13472   isMappedNone    NA      NA      NA      NA
7734    isMappedwt      mqwtpass        NA      NA      NA
152     isMappeddel     NA      mqdelpass       NA      NA
"""

# 
# can't identify the origin of these reads:  822423
# 
# for those that can be placed uniquely:
#   wt: 689502 + 32436 + 7734 =  729672  (53.7%)
#   del: 595223 + 33068 + 152 = 628443   (46.3%)
#

