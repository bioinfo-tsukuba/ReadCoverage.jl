
"""
    load_transcript(path_bed12::String)

Loads transcripts information from a BED12-format file

# Arguments
- `path_bed12::String`: the path of the file
"""
function load_transcript(path_bed12::String)
	transcripts = Array{BED.Record,1}()
	open(BED.Reader, path_bed12) do bed12_reader
		for record::BED.Record in bed12_reader
			push!(transcripts, record)
		end
	end
	return transcripts
end

