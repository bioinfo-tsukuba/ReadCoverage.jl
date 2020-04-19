export absolute_genebodycoverage

"""
absolute_genebodycoverage(path_bam::String, path_bed12::String, output_prefix::String; bin_size::Int=100)
Todo.
Arguments
---------
- `path_bam`: Todo
- `path_bed12`: Todo
- `output_prefix`: Todo
- `bin_size`: Todo
"""
function absolute_genebodycoverage(path_bam::String, path_bed12::String, output_prefix::String; bin_size::Int=100)
	path_out = output_prefix * ".absoluteGeneBodyCoverage.txt"
	println(@sprintf "bam: %s\nbed12: %s\npath_out: %s\nbin_size: %d" path_bam path_bed12 path_out bin_size)

	# File check
	if !isfile(path_bam)
		error(@sprintf "No such file: %s\n" path_bam)
	end
	if !isfile(path_bed12)
		error(@sprintf "No such file: %s\n" path_bed12)
	end
	if !uv_access_writable(dirname(output_prefix))
		error(@sprintf "Output files are not writable to: %s\n" (dirname(output_prefix) != "" ? dirname(output_prefix) : "."))
	end

	# Load transcripts information (BED12)
	transcripts = Array{BED.Record,1}()
	open(BED.Reader, path_bed12) do bed12_reader
		for record::BED.Record in bed12_reader
			push!(transcripts, record)
		end
	end
	println(@sprintf "transcripts: %d" length(transcripts))

	# Calculate transcript length
	transcripts_length = zeros(Int, length(transcripts))
	for i in 1:length(transcripts)
		transcripts_length[i] = sum(BED.blocksizes(transcripts[i]))
	end

	# Define L_max (Maximum length of transcript)
	L_max = 0
	L_max = maximum(transcripts_length)
	println(@sprintf("L_max: %d", L_max))

	# Calculate K
	K = Int(ceil(L_max/bin_size))
	println(@sprintf("K: %d", K))

	# Calculate number of transcripts in each bin
	N_transcript = zeros(Int, K)
	for L in transcripts_length
		N_transcript[Int(ceil(L/bin_size))] += 1
	end
	for k in (K-1):-1:1
		N_transcript[k] += N_transcript[k+1]
	end

	# Generate a L_max-length array for read coverage
	cov = zeros(L_max)

	# Generate a K-length array for absolute coverage
	abcov = zeros(K)

	# Open BAM file
	open(BAM.Reader, path_bam, index=path_bam*".bai") do bam_reader
		# For each transcript
		println("Calculate absolute read coverage for each transcript...")

		for t in transcripts
			# Get blockSizes (= exon lengths) and blockStarts (=start position of exons relative to chromStart of t)
			blockSizes = BED.blocksizes(t)
			L = sum(blockSizes)
			blockStarts = BED.blockstarts(t)

			# Strand
			BED.strand(t)


			# Initialize cov
			for p in 1:length(cov)
				cov[p] = 0
			end

			# Get the read coverage on each exon and save it to `cov`
			# Note: Read coverage is saved in `cov` in a 'Left justified' manner,
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

				cov[pos_cov] = readcoverage_bam_base(bam_reader, BED.chrom(t), exon_start, exon_end)
			end

			# Calculate count per bin
			# (Distance from 3'-end is binned.)
			# (k-th bin: [(k-1)*b+1, k*b])
			for k in 1:K # for each bin
				if k*bin_size > L # Case when k-th bin exceeds boundary of transcript
					if sum(cov[1:(L-((k-1)*bin_size+1)+1)]) > 0
						abcov[k] += 1
						break
					end
				else
					if sum(cov[(L-k*bin_size+1):(L-((k-1)*bin_size+1)+1)]) > 0
						abcov[k] += 1
					end
				end
			end
		end
	end

	# Divide count by the number of transcripts for each bin
	for k in 1:K
		abcov[k] = abcov[k] / N_transcript[k]
	end

	# Write output (TSV) (Sample_ID, Bin, Coverage)
	println("Write results...")
	binNumbers = collect(1:K)
	binStarts = (binNumbers .- 1) .* bin_size .+ 1
	binEnds = binNumbers .* bin_size
	open(path_out, "w") do io
           writedlm(io, [binNumbers binStarts binEnds abcov N_transcript], '\t')
	end
	
	# Save plot
	plot_absolute_coverage(abcov, out_path=output_prefix * ".absoluteGeneBodyCoverage.pdf")

    return(abcov)
end
