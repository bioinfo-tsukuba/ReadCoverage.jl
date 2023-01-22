module ReadCoverage

using GenomicFeatures
using XAM
using BED
using BioAlignments
using DelimitedFiles
using Printf
using Plots
using Statistics
using GZip

include("relative_genebodycoverage.jl")
include("readcovearge_bam.jl")
include("absolute_genebodycoverage.jl")
include("load_transcript.jl")
include("plot.jl")
include("utils.jl")
include("vdoc.jl")

end # module