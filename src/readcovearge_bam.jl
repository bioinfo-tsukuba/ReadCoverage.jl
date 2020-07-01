export readcoverage_bam
export readcoverage_transcript_bam

"""
	readcoverage_bam(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64;
					output_prefix::String = "")

Calculates a read coverage for a genomic interval.

# Arguments
---------
- `path_bam::String`: the path to a BAM file.
- `chrom::String`: the chromosome name for a genomic interval.
- `leftpos::Int64`: the left position for a genomic interval.
- `rightpos::Int64`: the right position for a genomic interval.
- `output_prefix::String`: the prefix for output files. If this keyword is not specified or set to "" (defalut), no output files are saved.
"""
function readcoverage_bam(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64; output_prefix::String = "")
	if output_prefix != ""
		path_out = output_prefix * ".readCoverage.txt"
		path_out_plot = output_prefix * ".readCoverage.pdf"
		println(@sprintf "Start calculating read coverage...\n- bam: %s\n- output: %s\n          %s\n" path_bam path_out path_out_plot)
	end

	# File check
	if !isfile(path_bam)
		error(@sprintf "No such file: %s\n" path_bam)
	end
	if output_prefix != ""
		if !uv_access_writable(dirname(output_prefix))
			error(@sprintf "Output files are not writable to: %s\n" (dirname(output_prefix) != "" ? dirname(output_prefix) : "."))
		end
	end

	# Define the output read coverage container
	coverage = zeros(Int, rightpos - leftpos + 1)

	# Caclutate the read coverage
	open(BAM.Reader, path_bam, index=path_bam*".bai") do bam_reader

		# Check wheteher the region exists in BAM
		if ! exists_region_bam(bam_reader, chrom, leftpos, rightpos)
			error("The region does not exist in the BAM file.")
		end

        coverage = readcoverage_bam_base(bam_reader, chrom, leftpos, rightpos)
	end
	
	if output_prefix != ""

		# Write output (TSV) (Position, Coverage)
		println("Write results...")
		bin_numbers = collect(leftpos:rightpos)
		open(path_out, "w") do io
			writedlm(io, [bin_numbers coverage], '\t')
		end
		
		# Save plot
		plot_read_coverage(coverage, out_path = path_out_plot)

		# Message
		println(@sprintf "Finished! Check output files:\n- %s\n- %s" path_out path_out_plot)
	end

	return(coverage)
end


"""
	readcoverage_bam_base(bam_reader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)

Generates the read coverage for a genomic interval.
Both `leftpos` and `rightpos` coordinates are assumed to be 1-based.
"""
function readcoverage_bam_base(bam_reader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)
	# Define output read coverage container
	coverage = zeros(Int, rightpos - leftpos + 1)

	# Evaluate each BAM record
	for record in eachoverlap(bam_reader, chrom, leftpos:rightpos)

        # Skip a BAM record for an unmapped read (if pileupread.alignment.is_unmapped:continue)
		if ! BAM.ismapped(record)
			continue
        end

        # Compatibility with RSeQC genebody_coverage.py
        ## if pileupread.alignment.is_secondary:continue
        if BAM.flag(record) & SAM.FLAG_SECONDARY != 0
            continue
        end

        ## if pileupread.alignment.is_qcfail:continue 
        if BAM.flag(record) & SAM.FLAG_QCFAIL != 0
            continue
        end

        ## if pileupread.alignment.is_duplicate:continue
        if BAM.flag(record) & SAM.FLAG_DUP != 0
            continue
        end           
        
		read_leftpos = BAM.position(record)
		cigarRle = BAM.cigar_rle(record)
		offset = 0

		# Decode CIGAR string to get coverage of split-aligned reads
		for i in 1:length(cigarRle[1])
			if cigarRle[1][i] == OP_MATCH
				if leftpos <= read_leftpos + offset <= rightpos
					for j in 1:min(cigarRle[2][i], rightpos - read_leftpos - offset + 1)
						coverage[read_leftpos + offset + j - leftpos] += 1
					end
				elseif read_leftpos + offset < leftpos
					for j in (leftpos-(read_leftpos + offset)+1):min(cigarRle[2][i], rightpos - read_leftpos - offset + 1)
						coverage[read_leftpos + offset + j - leftpos] += 1
					end
				end
			end
			offset += cigarRle[2][i]
		end
	end

	return(coverage)
end


"""
	readcoverage_transcript_bam(bam_reader::BAM.Reader, t::BED.Record)

Returns read coverage for a transcript.
"""
function readcoverage_transcript_bam(bam_reader::BAM.Reader, t::BED.Record)

	# Check if the region of the transcript exists in BAM
	if ! exists_region_bam(bam_reader, BED.chrom(t), BED.chromstart(t), BED.chromend(t))
		error("The transcript region does not exist in the BAM file.")
	end

	# Get block_sizes (= exon lengths) and block_starts (=start position of exons relative to chromStart of t)
	block_sizes = BED.blocksizes(t)
	transcript_length = sum(block_sizes)
	block_starts = BED.blockstarts(t)

	coverage = zeros(Float64, transcript_length)

	# Get the read coverage on each exon and save it to `cov`
	# Note: Read coverage is saved in `cov` in a 'Left justified' manner,
	#       and left is 5'-end.
	# If the strand of a transcript is '-' STRAND_NEG, read coverage was reversed.
	p_start = 1
	p_end = 0
	s_g = BED.chromstart(t) + 1
	e_g = 0

	for i in 1:length(block_sizes) # For each exon
	
		# Calculate start and end positions of an exon on the genome
		s_g = block_starts[i] + s_g
		e_g = s_g + block_sizes[i] - 1

		# Calculate start and end positions of an exon on mRNA
		p_start = p_end + 1
		p_end = p_start + block_sizes[i] - 1

		# Get read coverage and save it
		if BED.strand(t) == STRAND_POS
			pos_cov = p_start:p_end
		elseif BED.strand(t) == STRAND_NEG
			pos_cov = (transcript_length-p_start+1):-1:(transcript_length-p_end+1)
		else # STRAND_NA, STRAND_BOTH
			continue
		end

		coverage[pos_cov] = readcoverage_bam_base(bam_reader, BED.chrom(t), s_g, e_g)
	end

	return(coverage)
end


"""
For test
"""
function bamToCoverage_cigarAware_position(bam_reader::BAM.Reader, chrom::String, pos::Int64)
	coverage = 0

	# Evaluate each BAM record
	for record in eachoverlap(bam_reader, chrom, pos:pos)
		# Skip a BAM record for an unmapped read
		if ! BAM.ismapped(record)
			continue
		end

		read_leftpos = BAM.position(record)
		cigarRle = BAM.cigar_rle(record)
		offset = 0

		# Decode CIGAR string to get coverage of split-aligned reads
		for i in 1:length(cigarRle[1])
			if cigarRle[1][i] == OP_MATCH
				if read_leftpos + offset <= pos && pos <= read_leftpos + offset + cigarRle[2][i] - 1
					coverage += 1
				end
			end
			offset += cigarRle[2][i]
		end
	end

	return(coverage)
end


"""
	exists_region_bam(bam_reader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)

Returns true if a region chrom:leftpos-rightpos exists in a BAM file for which a reader is opened.
"""
function exists_region_bam(bam_reader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)::Bool
	if in(chrom, bam_reader.refseqnames) &&
			leftpos <= rightpos &&
			leftpos >= 0 &&
			rightpos <= bam_reader.refseqlens[findfirst(isequal(chrom), bam_reader.refseqnames)]
		return true
	else
		return false
	end
end
