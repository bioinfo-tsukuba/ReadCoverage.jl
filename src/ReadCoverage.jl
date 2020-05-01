module ReadCoverage

using GenomicFeatures
using XAM
using BioAlignments
using DelimitedFiles
using Printf
using Plots
using PyPlot

include("relative_genebodycoverage.jl")
include("readcovearge_bam.jl")
include("absolute_genebodycoverage.jl")
include("load_transcript.jl")
include("plot.jl")
include("utils.jl")

end # module