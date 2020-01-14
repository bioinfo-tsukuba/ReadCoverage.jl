

# Ref: GenomicFeatures.jl https://biojulia.net/GenomicFeatures.jl/latest/io/bed/
# Ref strand!! https://github.com/BioJulia/GenomicFeatures.jl/blob/8fc34ff680f5e742e25d2bd3d4722cb33fbe3cd5/src/strand.jl
# Future TODO: Test
# Future TODO: How to make subset of 'expressed' transcripts (Give list of transcript id?)
# Ref: Plot https://qiita.com/yuifu/items/b2031d0367da878fddae#bioalignments
# Future TODO: Consider Stranded RNA-seq

using GenomicFeatures
using BioAlignments
using DelimitedFiles
using Printf


# CAUTION: No strandness
function bamToCoverage_cigarAware_base(bamReader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)
	coverage = zeros(Int, rightpos - leftpos + 1)
	for record in eachoverlap(bamReader, chrom, leftpos:rightpos)
		# skip
		if ! BAM.ismapped(record)
			continue
		end

		readLeftPos = BAM.position(record)
		cigarRle = BAM.cigar_rle(record)
		offset = 0

		for i in 1:length(cigarRle[1])
			if cigarRle[1][i] == OP_MATCH
				if leftpos <= readLeftPos + offset <= rightpos
					for j in 1:min(cigarRle[2][i], rightpos - readLeftPos - offset + 1)
						coverage[readLeftPos + offset + j - leftpos] += 1
					end
				elseif readLeftPos + offset < leftpos
					for j in (leftpos-(readLeftPos + offset)+1):min(cigarRle[2][i], rightpos - readLeftPos - offset + 1)
						coverage[readLeftPos + offset + j - leftpos] += 1
					end
				end
			end
			offset += cigarRle[2][i]
		end
	end
	return(coverage)
end

function bamToCoverage_cigarAware(pathBam::String, chrom::String, leftpos::Int64, rightpos::Int64)
    open(BAM.Reader, pathBam, index=pathBam*".bai") do reader
        return(bamToCoverage_cigarAware_base(reader, chrom, leftpos, rightpos))
    end
end

"""
	abCov(path_bam::String, path_bed12::String, path_out::String; bin_size::Int=100)
Todo.
Arguments
---------
- `path_bam`: Todo
- `path_bed12`: Todo
- `path_out`: Todo
- `bin_size`: Todo
"""
function abCov(path_bam::String, path_bed12::String, path_out::String; bin_size::Int=100)

	println(@sprintf "bam: %s\nbed12: %s\npath_out: %s\nbin_size: %d" path_bam path_bed12 path_out bin_size)

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
			s_g = BED.chromstart(t) + 1
			e_g = 0
			for i in 1:length(blockSizes) # For each exon
				# Calculate start and end positions of an exon on the genome
				s_g = blockStarts[i] + s_g
				e_g = s_g + blockSizes[i] - 1

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

				cov[pos_cov] = bamToCoverage_cigarAware_base(bam_reader, BED.chrom(t), s_g, e_g)
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

    return(abcov)
end
