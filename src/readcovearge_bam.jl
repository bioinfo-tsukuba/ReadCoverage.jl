"""
	Generates read coverage for a genomic interval.
	Both `leftpos` and `rightpos` coordinates are assumed to be 1-based.
"""
function bamToCoverage_cigarAware_base(bamReader::BAM.Reader, chrom::String, leftpos::Int64, rightpos::Int64)
	# Define output read coverage container
	coverage = zeros(Int, rightpos - leftpos + 1)

	# Evaluate each BAM record
	for record in eachoverlap(bamReader, chrom, leftpos:rightpos)
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
        
		readLeftPos = BAM.position(record)
		cigarRle = BAM.cigar_rle(record)
		offset = 0

		# Decode CIGAR string to get coverage of split-aligned reads
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


"""
@deprectate
"""
function bamToCoverage_cigarAware(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64)
    open(BAM.Reader, path_bam, index=path_bam*".bai") do reader
        return(bamToCoverage_cigarAware_base(reader, chrom, leftpos, rightpos))
    end
end





"""
	Returns read coverage for a transcript
"""
function coverage_transcript(bam_reader::BAM.Reader, t::BED.Record)
	# Get blockSizes (= exon lengths) and blockStarts (=start position of exons relative to chromStart of t)
	blockSizes = BED.blocksizes(t)
	L = sum(blockSizes)
	blockStarts = BED.blockstarts(t)

	coverage = zeros(Float64, L)

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

		coverage[pos_cov] = bamToCoverage_cigarAware_base(bam_reader, BED.chrom(t), s_g, e_g)
	end

	return(coverage)
end




"""
@deprectate
Test
"""
function bamToCoverage_cigarAware_position(bamReader::BAM.Reader, chrom::String, pos::Int64)
	coverage = 0

	# Evaluate each BAM record
	for record in eachoverlap(bamReader, chrom, pos:pos)
		# Skip a BAM record for an unmapped read
		if ! BAM.ismapped(record)
			continue
		end

		readLeftPos = BAM.position(record)
		cigarRle = BAM.cigar_rle(record)
		offset = 0

		# Decode CIGAR string to get coverage of split-aligned reads
		for i in 1:length(cigarRle[1])
			if cigarRle[1][i] == OP_MATCH
				if readLeftPos + offset <= pos && pos <= readLeftPos + offset + cigarRle[2][i] - 1
					coverage += 1
				end
			end
			offset += cigarRle[2][i]
		end
	end

	return(coverage)
end


