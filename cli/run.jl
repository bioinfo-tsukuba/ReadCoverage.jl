using ArgParse
using ReadCoverage

println("####################\n   ReadCoverage.jl\n####################")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "relcov"
            help = "Calculates relative gene body coverage from a BAM file.\n(Almost same as genebody_coverage.py in RSeQC.)"
            action = :command        
        "abcov"
            help = "Cacluates absolute gene body coverage from a BAM file."
            action = :command
        "coverage"    
            help = "Cacluates absolute gene body coverage from a BAM file."
            action = :command
    end

    @add_arg_table! s["relcov"] begin
        "path_bam"
            help = "Path to a BAM file."
            arg_type = String
            required = true
        "path_bed12"
            help = "Path to a reference gene model file with BED12 format."
            arg_type = String
            required = true
        "output_prefix"
            help = "Prefix for output files. If not specified or set to \"\" (defalut), no output files are saved."
            arg_type = String
            required = true
        "--max_depth"
            help = "Maximum depth for a position. Depth more than `max_depth` is cut to `max_depth`. If this keyword is set to 0 (defalut), no cut is occurred."
            arg_type = Int
            default = 0   
    end

    @add_arg_table! s["abcov"] begin
        "path_bam"
            help = "Path to a BAM file."
            arg_type = String
            required = true
        "path_bed12"
            help = "Path to a reference gene model file with BED12 format."
            arg_type = String
            required = true
        "output_prefix"
            help = "refix for output files. If not specified or set to \"\" (defalut), no output files are saved."
            arg_type = String
            required = true
        "--bin_size"
            help = "Bin size for transcripts (Defalut: 100)."
            arg_type = Int
            default = 100 
    end
    @add_arg_table! s["coverage"] begin
        "path_bam"
            help = "Path to a BAM file."
            arg_type = String
            required = true
        "chrom"
            help = "Chromosome name for a genomic interval."
            arg_type = String
            required = true
        "leftpos"
            help = "Left position for a genomic interval."
            arg_type = Int64
            required = true
        "rightpos"
            help = "Right position for a genomic interval."
            arg_type = Int64
            required = true
        "output_prefix"
            help = "Prefix for output files. If not specified or set to \"\" (defalut), no output files are saved."
            arg_type = String
            required = true
    end
    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()

    if parsed_args["%COMMAND%"] == "relcov"
        relative_genebodycoverage(
            parsed_args["relcov"]["path_bam"], 
            parsed_args["relcov"]["path_bed12"], 
            output_prefix = parsed_args["relcov"]["output_prefix"],
            max_depth=parsed_args["relcov"]["max_depth"]
        );
    elseif parsed_args["%COMMAND%"] == "abcov"
        absolute_genebodycoverage(
            parsed_args["abcov"]["path_bam"], 
            parsed_args["abcov"]["path_bed12"], 
            output_prefix = parsed_args["abcov"]["output_prefix"],
            bin_size=parsed_args["abcov"]["bin_size"]
        );
    elseif parsed_args["%COMMAND%"] == "coverage"
        readcoverage_bam(
            parsed_args["coverage"]["path_bam"],
            parsed_args["coverage"]["chrom"],
            parsed_args["coverage"]["leftpos"],
            parsed_args["coverage"]["rightpos"],
            output_prefix = parsed_args["coverage"]["output_prefix"]
        );
    end

end

main()
