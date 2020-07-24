export vdoc

"""
    vdoc(array_path_bam::Array{String,1}, path_bed12::String;
					output_prefix::String = "")

Calculates Variability of Depth of Coverage (VDoC) score across samples.

A VDoC score is defined for a transcript as follows:
1. For each base on a transcript, a depth of coverage is calculated. 
2. For each base, the coefficient of variation (CV) of the depths of coverage is calculated
   across samples.
3. The average of CVs across the transcript is calculated, which we designate as a VDoC score.

# Arguments
---------
- `array_path_bam::Array{String,1}`: Array of paths to BAM files.
- `path_bed12::String`: Path to a reference gene model file with BED12 format.
- `output_prefix::String`: Prefix for output files. If this keyword is not specified or set to
  "" (defalut), no output files are saved.
"""
function vdoc(array_path_bam::Array{String,1}, path_bed12::String;
					output_prefix::String = "")

	if output_prefix != ""
		path_out = output_prefix * ".vdoc.txt"
		path_out_plot = output_prefix * ".vdoc.pdf"
		println(@sprintf "Start calculating absolute gene body coverage...\n- bam: %s\n- bed12: %s\n- output: %s\n          %s\n" join(array_path_bam, ", ") path_bed12 path_out path_out_plot)
	end

	# File check
	for i in 1:length(array_path_bam)
		path_bam = array_path_bam[i]
		if !isfile(path_bam)
			error(@sprintf "No such file: %s\n" path_bam)
		end
	end
	if !isfile(path_bed12)
		error(@sprintf "No such file: %s\n" path_bed12)
	end
	if output_prefix != ""
		if !uv_access_writable(dirname(output_prefix))
			error(@sprintf "Output files are not writable to: %s\n" (dirname(output_prefix) != "" ? dirname(output_prefix) : "."))
		end
	end

	# Load transcripts								
	transcripts = load_transcript(path_bed12)
	println(@sprintf "transcripts: %d" length(transcripts))

	# Define containers of results
	array_vdoc = zeros(length(transcripts))
	array_transcript_length = zeros(Int64, length(transcripts))
	array_transcript_name = Array{String,1}(undef, length(transcripts))

	# Caclutate the read coverage
	for (j, t) in enumerate(transcripts)
		# Define matrix
		transcript_length = sum(BED.blocksizes(t))
		coverages = zeros(Int, length(array_path_bam), transcript_length)

		for (i, path_bam) in enumerate(array_path_bam)
			open(BAM.Reader, path_bam, index=path_bam*".bai") do bam_reader
				# Load coverage
				coverages[i,:] = readcoverage_transcript_bam(bam_reader::BAM.Reader, t::BED.Record)
			end
		end
		# Calculate mean of CV
		mean_coverages = mean(coverages, dims=1)
		std_coverages = std(coverages, dims=1)
		cv_coverages = std_coverages ./ mean_coverages
		vdoc = mean(cv_coverages[mean_coverages .> 0])

		array_vdoc[j] = vdoc
		array_transcript_length[j] = transcript_length
		array_transcript_name[j] = BED.name(t)
	end

	if output_prefix != ""
		# Write output (TSV) (Sample_ID, Bin, Coverage)
		println("Write results...")
		open(path_out, "w") do io
			writedlm(io, [array_transcript_name array_transcript_length array_vdoc], '\t')
		end

		# Message
		println(@sprintf "Finished! Check output files:\n- %s" path_out)
	end
	
    return(array_transcript_name, array_transcript_length, array_vdoc)
end
