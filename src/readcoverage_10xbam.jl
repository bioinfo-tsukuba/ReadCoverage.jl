export read_10xbam
export plot_10xbam


"""
"""
function read_10xbam(path_bam::String, chrom::String, leftpos::Int64, rightpos::Int64, path_cellbarcodes::String)

    # Load cell barcode information
    cellbarcodes = read_cellbarcodes(path_cellbarcodes)
    n_barcodes = length(keys(cellbarcodes))

    # Make output
    outmaterix = zeros(Int64, n_barcodes, rightpos - leftpos + 1)

    # Open a BAM file
    bam_reader = open(BAM.Reader, path_bam, index = path_bam * ".bai")

    cell_index = 0

    # Iterate over records overlapping the specified region
	#try
		# Evaluate each BAM record
        for record in eachoverlap(bam_reader, chrom, leftpos:rightpos)
			# Skip a BAM record for an unmapped read (if pileupread.alignment.is_unmapped:continue)
			if ! BAM.ismapped(record)
				continue
			end

			# # Compatibility with RSeQC genebody_coverage.py
			# ## if pileupread.alignment.is_secondary:continue
			# if BAM.flag(record) & SAM.FLAG_SECONDARY != 0
			# 	continue
			# end

			# ## if pileupread.alignment.is_qcfail:continue 
			# if BAM.flag(record) & SAM.FLAG_QCFAIL != 0
			# 	continue
			# end

			# ## if pileupread.alignment.is_duplicate:continue
			# if BAM.flag(record) & SAM.FLAG_DUP != 0
			# 	continue
			# end
			
			read_leftpos = BAM.position(record)
			cigarRle = BAM.cigar_rle(record)
            offset = 0           

            # If the cell barcode is defined
            try
                if haskey(cellbarcodes, BAM.auxdata(record)["CB"])
                    cell_index = cellbarcodes[BAM.auxdata(record)["CB"]]
                else
                    continue
                end
            catch e
                continue
            end
            # println(BAM.refname(record), ':', BAM.position(record))
            # return(record)

			# Decode CIGAR string to get coverage of split-aligned reads
			for i in 1:length(cigarRle[1])
				if cigarRle[1][i] == OP_MATCH
					if leftpos <= read_leftpos + offset <= rightpos
						for j in 1:min(cigarRle[2][i], rightpos - read_leftpos - offset + 1)
							outmaterix[cell_index, read_leftpos + offset + j - leftpos] += 1
						end
					elseif read_leftpos + offset < leftpos
						for j in (leftpos-(read_leftpos + offset)+1):min(cigarRle[2][i], rightpos - read_leftpos - offset + 1)
							outmaterix[cell_index, read_leftpos + offset + j - leftpos] += 1
						end
					end
				end
				offset += cigarRle[2][i]
			end
		end
	# catch e
	#	println("Catched and ignored an error: ", typeof(e))
    #end
    
    # Close the BAM file.
    close(bam_reader)

    return(outmaterix)
end


"""
"""
function read_cellbarcodes(path_cellbarcode::String)::Dict{String,Int64}
    cellbarcodes = Dict{String,Int64}()
    GZip.open(path_cellbarcode, "r") do f
        for (i, barcode) in enumerate(readlines(f))
            cellbarcodes[barcode] = i
        end
    end
    return(cellbarcodes)
end

"""
"""
function hastag(record::XAM.BAM.Record, tag::AbstractString)
    for aux in XAM.BAM.auxdata(record)
        if aux[1] == tag
            return true
        end
    end
    return false
end
