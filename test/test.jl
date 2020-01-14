include("../src/AbCov.jl")


path_bam = "/Users/haruka/RamDA_10pg_12.4M.1.uniq.proper_pair.bam"
path_bed12 = "/Users/haruka/Dropbox (Personal)/research/2019/10_abCov/AbCov/test/test_data/gencode.vM15.primary_assembly.annotation.protein_coding.head.bed"
path_out = "/Users/haruka/Dropbox (Personal)/research/2019/10_abCov/AbCov/test/test_data/out.txt"

abcov = abCov(path_bam, path_bed12, path_out; bin_size=100)


using Plots
x = reverse(collect(1:length(abcov)) .* 100)
Plots.plot(x, abcov)


using PyPlot
x = reverse(collect(1:length(abcov)) .* 100)
PyPlot.plot(x, abcov, linewidth=2.0)
gcf()
clf()


############
chrom = "chr19"
leftpos = 3205815
rightpos = 3206000
cov = bamToCoverage_cigarAware(path_bam, chrom, leftpos, rightpos)
Plots.plot(leftpos:rightpos, cov)
