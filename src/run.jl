include("AbCov.jl")


path_bam = ARGS[1]
path_bed12 = ARGS[2]
output_prefix = ARGS[3]

max_depth = 0
if length(ARGS) > 3
    max_depth = parse(Int, ARGS[4])    
end

coverage = relative_genebodycoverage(path_bam, path_bed12, output_prefix, max_depth=max_depth);
