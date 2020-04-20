using ReadCoverage

path_bam = "examples/RamDA_72h_A09.uniq.q40.bam"
path_bed12 = "examples/gencode.vM15.primary_assembly.annotation.protein_coding.head.bed"
out_prefix = "examples/out/RamDA_72h_A09.uniq.q40"


############
chrom = "chr19"
leftpos = 3205000
rightpos = 3207000
cov = readcoverage_bam(path_bam, chrom, leftpos, rightpos)
cov = readcoverage_bam(path_bam, chrom, leftpos, rightpos, output_prefix=out_prefix)
plot_read_coverage(cov)


################################
# Absolute genebody coverage
abcov = absolute_genebodycoverage(path_bam, path_bed12)
abcov = absolute_genebodycoverage(path_bam, path_bed12, output_prefix=out_prefix, bin_size=100)
plot_absolute_coverage(abcov)


################################
# Absolute genebody coverage
relcov = relative_genebodycoverage(path_bam, path_bed12);

max_depth = 0
relcov = relative_genebodycoverage(path_bam, path_bed12, output_prefix=out_prefix, max_depth=max_depth);
plot_relative_coverage(relcov)
