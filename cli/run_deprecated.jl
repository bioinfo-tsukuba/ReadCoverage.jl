using Fire
using ReadCoverage

"""
    relcov(path_bam::String, path_bed12::String, out::String; max_depth::Int=0)

Calculates relative gene body coverage from a BAM file.
(Almost same as genebody_coverage.py in RSeQC.)
"""
@main function relcov(path_bam::String, path_bed12::String, out::String; max_depth::Int=0)
    relative_genebodycoverage(path_bam, path_bed12, output_prefix=out, max_depth=max_depth);
end

"""
    abcov(path_bam::String, path_bed12::String, out::String; bin_size::Int=100)

Cacluates absolute gene body coverage from a BAM file
"""
@main function abcov(path_bam::String, path_bed12::String, out::String; bin_size::Int=100)
    absolute_genebodycoverage(path_bam, path_bed12, output_prefix=out, bin_size=bin_size);
end

"""
    coverage(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64, out::String)
    
Cacluates absolute gene body coverage from a BAM file
"""
@main function coverage(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64, out::String)
    readcoverage_bam(path_bam, chrom, leftpos, rightpos, output_prefix=out);
end
