export absolute_genebodycoverage

"""
    absolute_genebodycoverage(path_bam::String, path_bed12::String;
							output_prefix::String = "", bin_size::Int=100)

Calculates absolute gene body coverage.

# Arguments
---------
- `path_bam::String`: Path to a BAM file.
- `path_bed12::String`: Path to a reference gene model file with BED12 format.
- `output_prefix::String`: Prefix for output files. If this keyword is not specified or set to "" (defalut), no output files are saved.
- `bin_size::Int`: Bin size for transcripts (Defalut: 100).
"""
function absolute_genebodycoverage(path_bam::String, path_bed12::String;
									output_prefix::String = "", bin_size::Int=100)
	if output_prefix != ""
		path_out = output_prefix * ".absoluteGeneBodyCoverage.txt"
		path_out_plot = output_prefix * ".absoluteGeneBodyCoverage.pdf"
		println(@sprintf "Start calculating absolute gene body coverage...\n- bam: %s\n- bed12: %s\n- output: %s\n          %s\n- bin_size: %d\n" path_bam path_bed12 path_out path_out_plot bin_size)
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

	# Load transcripts information (BED12)
	transcripts = Array{BED.Record,1}()
	open(BED.Reader, path_bed12) do bed12_reader
		for record::BED.Record in bed12_reader
			push!(transcripts, record)
		end
	end
	println(@sprintf "transcripts: %d" length(transcripts))

	# Calculate L_max (Maximum length of transcript)
	L_max = maximum(sum(BED.blocksizes(t)) for t in transcripts)
	println(@sprintf("L_max: %d", L_max))

	# Calculate K
	K = Int(ceil(L_max/bin_size))
	println(@sprintf("K: %d", K))

	# Calculate the number of transcripts in each bin
	N_transcript = zeros(Int, K)

	# Generate a L_max-length array for read coverage
	coverage = zeros(L_max)

	# Generate a K-length array for absolute coverage
	abcov = zeros(K)

	# Open BAM file
	open(BAM.Reader, path_bam, index=path_bam*".bai") do bam_reader
	
		# For each transcript
		println("Calculate absolute read coverage for each transcript...")

		for t in transcripts
			# Skip the transcript whose region does not exist in BAM
			if ! exists_region_bam(bam_reader, BED.chrom(t), BED.chromstart(t), BED.chromend(t))
				continue
			end

			# Get blockSizes (= exon lengths) and blockStarts (=start position of exons relative to chromStart of t)
			blockSizes = BED.blocksizes(t)
			L = sum(blockSizes)
			blockStarts = BED.blockstarts(t)

			# Count the number of transcripts whose 5'-end in each bin
			N_transcript[Int(ceil(sum(blockSizes)/bin_size))] += 1

			# Strand
			BED.strand(t)


			# Initialize coverage
			for p in 1:length(coverage)
				coverage[p] = 0
			end

			# Get the read coverage on each exon and save it to `coverage`
			# Note: Read coverage is saved in `coverage` in a 'Left justified' manner,
			#       and left is 5'-end.
			# If the strand of a transcript is '-' STRAND_NEG, read coverage was reversed.
			p_start = 1
			p_end = 0
			exon_start = BED.chromstart(t) # Returns 1-based coordinate(though 0-based in BED12 format)
			exon_end = 0
			for i in 1:length(blockSizes) # For each exon
				# Calculate start and end positions of an exon on the genome
				exon_start = blockStarts[i] - 1 + exon_start  # Returns 1-based coordinate(though 0-based in BED12 format)
				exon_end = exon_start + blockSizes[i] - 1

				# Calculate start and end positions of an exon on mRNA
				p_start = p_end + 1
				p_end = p_start + blockSizes[i] - 1

				# Get read coverage and save it
				if BED.strand(t) == STRAND_POS
					pos_cov = p_start:p_end
				elseif BED.strand(t) == STRAND_NEG
					pos_cov = (L-p_start+1):-1:(L-p_end+1)
				else # STRAND_NA, STRAND_BOTH
					continue
				end

				coverage[pos_cov] = readcoverage_bam_base(bam_reader, BED.chrom(t), exon_start, exon_end)
			end

			# Calculate count per bin
			# (Distance from 3'-end is binned.)
			# (k-th bin: [(k-1)*b+1, k*b])
			for k in 1:K # for each bin
				if k*bin_size > L # Case when k-th bin exceeds boundary of transcript
					if sum(coverage[1:(L-((k-1)*bin_size+1)+1)]) > 0
						abcov[k] += 1
						break
					end
				else
					if sum(coverage[(L-k*bin_size+1):(L-((k-1)*bin_size+1)+1)]) > 0
						abcov[k] += 1
					end
				end
			end
		end
	end

	# Calculate number of transcripts in each bin
	for k in (K-1):-1:1
		N_transcript[k] += N_transcript[k+1]
	end

	# Divide count by the number of transcripts for each bin
	for k in 1:K
		abcov[k] = abcov[k] / N_transcript[k]
	end

	if output_prefix != ""
		# Write output (TSV) (Sample_ID, Bin, Coverage)
		println("Write results...")
		bin_numbers = collect(1:K)
		bin_starts = (bin_numbers .- 1) .* bin_size .+ 1
		bin_ends = bin_numbers .* bin_size
		open(path_out, "w") do io
			writedlm(io, [bin_numbers bin_starts bin_ends abcov N_transcript], '\t')
		end
		
		# Save plot
		plot_absolute_coverage(abcov, out_path=path_out_plot)

		# Message
		println(@sprintf "Finished! Check output files:\n- %s\n- %s" path_out path_out_plot)
	end
	
    return(abcov)
end
