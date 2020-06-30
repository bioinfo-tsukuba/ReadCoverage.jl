export relative_genebodycoverage

"""
    relative_genebodycoverage(path_bam::String, path_bed12::String; 
							output_prefix::String="", transcript_length_cut::Int=100, max_depth::Int=0, N_bin::Int = 100)

Calculates the relative gene body coverage.

# Arguments
---------
- `path_bam::String`: the path to a BAM file.
- `path_bed12::String`: the path to a reference gene model file with BED12 format.
- `output_prefix::String`: the prefix for output files. If this keyword is not specified or set to "" (defalut), no output files are saved.
- `transcript_length_cut`: the minimum length of transcripts to be used for calculation (Defalut: 100)
- `max_depth::Int`: the maximum depth for a position. Depth more than `max_depth` is cut to `max_depth`.
					If this keyword is set to 0 (defalut), no cut is occurred.
- `N_bin::Int`: the number of bins for a transcript (Defalut: 100).
"""
function relative_genebodycoverage(path_bam::String, path_bed12::String; 
									output_prefix::String="", transcript_length_cut::Int=100, max_depth::Int=0, N_bin::Int = 100)
	if output_prefix != ""
		path_out = output_prefix * ".geneBodyCoverage.txt"
		path_out_plot = output_prefix * ".geneBodyCoverage.pdf"
		println(@sprintf "Start calculating relative gene body coverage...\n- bam: %s\n- bed12: %s\n- output: %s\n          %s\n- transcript_length_cut: %d\n- max_depth: %d\n- N_bin: %d\n" path_bam path_bed12 path_out path_out_plot transcript_length_cut max_depth N_bin)
	end

	# File check
	if !isfile(path_bam)
		error(@sprintf "No such file: %s\n" path_bam)
	end
	if !isfile(path_bed12)
		error(@sprintf "No such file: %s\n" path_bed12)
	end
	if output_prefix != ""
		if !uv_access_writable(dirname(output_prefix))
			error(@sprintf "Output files are not writable to: %s\n" (dirname(output_prefix) != "" ? dirname(output_prefix) : "."))
		end
	end

	# Loads transcripts information from BED12-format file
	transcripts = load_transcript(path_bed12)
	println(@sprintf "transcripts: %d" length(transcripts))

	# Define relative gene body coverage output
	relcov = zeros(Float64, N_bin)

	# Open BAM file
	open(BAM.Reader, path_bam, index=path_bam*".bai") do bam_reader

		# For each transcript
		for t in transcripts

			# Skip 
			if sum(BED.blocksizes(t)) < transcript_length_cut
				continue
			end

			if BED.strand(t) != STRAND_POS && BED.strand(t) != STRAND_NEG  # STRAND_NA, STRAND_BOTH
				continue
			end

			# Skip the transcript whose region does not exist in BAM
			if ! exists_region_bam(bam_reader, BED.chrom(t), BED.chromstart(t), BED.chromend(t))
				continue
			end

			# Get read coverage of the percentile positions on a transcript
			if max_depth == 0
				relcov += coverage_transcript_percentile(bam_reader, t)
			else
				# The maximum number for (imcomplete) comaptibility with `--max-depth` parametr of `pileup()` used in RSeQC
				relcov += min.(max_depth, coverage_transcript_percentile(bam_reader, t))
			end
		end
    end

	if output_prefix != ""
		# Write result
		sample_id = replace(basename(path_bam), r"\.[^\.]+$" => s"")
		write_relative_genebodycoverage(relcov, path_out, sample_id)

		# Save plot
		plot_relative_coverage(relcov, out_path=path_out_plot)
	
		# Message
		println(@sprintf "Finished! Check output files:\n- %s\n- %s" path_out path_out_plot)
	end

	return relcov
end


"""
    coverage_transcript_percentile(bam_reader::BAM.Reader, t::BED.Record)

Gets read coverage on the percentile positions in a transcript
1 read for 1 transcript
"""
function coverage_transcript_percentile(bam_reader::BAM.Reader, t::BED.Record)

	# Get block_sizes  (= exon lengths) and block_starts (=start position of exons relative to chromStart of t)
	block_sizes  = BED.blocksizes(t)
	block_starts = BED.blockstarts(t) # Returns 1-based coordinate (though 0-based in BED12 format)
	transcript_start = BED.chromstart(t)  # Returns 1-based coordinate(though 0-based in BED12 format)
	transcript_end = BED.chromend(t)
	
	# Collect within-exon positions relative to the transcript start
	L = sum(block_sizes)
	exon_coordinates = zeros(Int, L)
	offset = 0
	for i in 1:length(block_sizes) # For each exon
		for j in 0:(block_sizes[i]-1)
			exon_coordinates[j+offset+1] = block_starts[i] + j
		end
		offset += block_sizes[i]
	end

	# Select percentile positions (resulting in 100 positions)
	percentile_positions = percentile_list(exon_coordinates)

	# Note: Read coverage is saved in `cov` in a 'Left justified' manner,
	#       and left is 5'-end.
	# If the strand of a transcript is '-' STRAND_NEG, read coverage was reversed.
	if BED.strand(t) == STRAND_NEG
		percentile_positions = reverse(percentile_positions)
	end
	
	return readcoverage_bam_base(bam_reader, BED.chrom(t), transcript_start, transcript_end)[percentile_positions]
end


"""
    write_relative_genebodycoverage(relative_coverage::Array{Float64,1}, path_out::String, sample_id::String)

Write result of `relative_genebodycoverage()`
"""
function  write_relative_genebodycoverage(relative_coverage::Array{Float64,1}, path_out::String, sample_id::String)
	open(path_out, "w") do fw
		# Header
		write(fw, "Percentile\t" * join(string.(collect(1:100)), "\t") * "\n")
		# Coverage
		write(fw, sample_id * "\t" * join(string.(relative_coverage), "\t") * "\n")
	end
end


"""
    convert2relativecoverage(converage::Array{Float64,1}, bin_number::Int)
Converts read coverage to relative gene body coverage

*future work to do (incomplete)
"""
function convert2relativecoverage(converage::Array{Float64,1}, bin_number::Int)
	relative_coverage = zeros(Float64, bin_number)
	L = length(coverage)


	for b in 1:bin_number
		
	end
end


"""
    percentile_list(N::Array{Int,1})

Finds the percentile of a list of values and returns the list of percentile of the values.

# Arguments
- 'N::Array{Int,1}': the list of values. Note N MUST BE already sorted.

Reference: https://github.com/MonashBioinformaticsPlatform/RSeQC/blob/cb42bd90afa8d131875b18d54f70821142b8ea79/rseqc/qcmodule/mystat.py#L156
"""
function percentile_list(N::Array{Int,1})
	if length(N) < 100
		return N
	end

	percentile = zeros(Int, 100)

	for i in 1:100
		k = (length(N)-1) * i/100.0
		f = floor(k)
		c = ceil(k)
		if f == c
			percentile[i] = N[Int(k)+1]
		else
			d0 = N[Int(f)+1] * (c-k)
			d1 = N[Int(c)+1] * (k-f)
			percentile[i] = Int(round(d0+d1))
		end
	end

	return percentile
end
