percentile_list(collect(1:200))

x = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450]
y = [103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 188, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247, 249, 276, 302, 304, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 363, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450]
percentile_list(x) == y

collect(1:100) == percentile_list(collect(0:100))





"""
Python

exon_starts = [100, 200, 300, 400]
exon_ends = [150, 250, 350, 450]

gene_all_base=[]
for st,end in zip(exon_starts,exon_ends):
    gene_all_base.extend(range(st+1,end+1))

print(gene_all_base)
>>> print(gene_all_base)
[101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450]

def percentile_list(N):
	"""
	Find the percentile of a list of values.
	@parameter N - is a list of values. Note N MUST BE already sorted.
    @return - the list of percentile of the values
	"""
	if not N:return None
	if len(N) <100: return N
	per_list=[]
	for i in range(1,101):
		k = (len(N)-1) * i/100.0
		f = math.floor(k)
		c = math.ceil(k)
		if f == c:
			per_list.append( int(N[int(k)])  )
		else:
			d0 = N[int(f)] * (c-k)
			d1 = N[int(c)] * (k-f)
			per_list.append(int(round(d0+d1)))	
    return per_list

import math

percentile_list (gene_all_base)
>>> percentile_list (gene_all_base)
[103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 188, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247, 249, 276, 302, 304, 306, 308, 310, 312, 314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336, 338, 340, 342, 344, 346, 348, 363, 402, 404, 406, 408, 410, 412, 414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 438, 440, 442, 444, 446, 448, 450]



>>> percentile_list([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100])
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
"""




############

path_bam="/data2/myoshimura/21_Nextflow_test/test2_githubcode_releaseTaskAfter_unstrandedSE/output_190527_nextflow_ramdaQC_SE_unstranded/170202_NS500723_0071_AH33LNBGX2/05_rseqc/gene_bodycoverage/NSRdT18_0127_C04_S28_R1_001_trim.sort.bam"
path_bed12="test/test_data/gencode.vM15.primary_assembly.annotation.protein_coding.head.bed"


transcript_length_cut = 100
@time transcripts = load_transcript(path_bed12);


reader = open(BAM.Reader, path_bam, index=path_bam*".bai")

# For each transcript
function unittest1(transcripts, reader, transcript_length_cut)
    println(length(transcripts))

    N_bin = 100
	# Define relative gene body coverage output
	relative_coverage = zeros(Float64, N_bin)
    for t in transcripts
        # 
        if sum(BED.blocksizes(t)) < transcript_length_cut
            continue
        end

        if BED.strand(t) != STRAND_POS && BED.strand(t) != STRAND_NEG  # STRAND_NA, STRAND_BOTH
            continue
        end

        # Get read coverage of the percentile positions on a transcript
        relative_coverage += coverage_transcript_percentile(reader, t)
    end
end

function unittest2(transcripts, reader, transcript_length_cut)
    println(length(transcripts))

    N_bin = 100
	# Define relative gene body coverage output
	relative_coverage = zeros(Float64, N_bin)
    for t in transcripts
        # 
        if sum(BED.blocksizes(t)) < transcript_length_cut
            continue
        end

        if BED.strand(t) != STRAND_POS && BED.strand(t) != STRAND_NEG  # STRAND_NA, STRAND_BOTH
            continue
        end

        # Get read coverage of the percentile positions on a transcript
        coverage_transcript_percentile(reader, t);
    end
end



@time unittest1(transcripts[1:2], reader, transcript_length_cut)

@time unittest1(transcripts[1:100], reader, transcript_length_cut);
100
 14.907310 seconds (538.96 k allocations: 49.787 MiB, 0.15% gc time)
# unittest1() の中にボトルネックがある

# relative_coverage への += は関係あるか? 少しは早くなるがあまり関係なさそう
@time unittest2(transcripts[1:100], reader, transcript_length_cut);
# 100
#  14.661962 seconds (555.47 k allocations: 50.516 MiB, 0.14% gc time)

# julia> @time unittest1(transcripts[1:100], reader, transcript_length_cut);
# 100
#   1.286662 seconds (49.66 k allocations: 9.811 MiB)


function unittest3(transcripts, reader, transcript_length_cut)
    println(length(transcripts))

    N_bin = 100
	# Define relative gene body coverage output
	relative_coverage = zeros(Float64, N_bin)
    for t in transcripts
        # 
        if sum(BED.blocksizes(t)) < transcript_length_cut
            continue
        end

        if BED.strand(t) != STRAND_POS && BED.strand(t) != STRAND_NEG  # STRAND_NA, STRAND_BOTH
            continue
        end

        blockSizes = BED.blocksizes(t)
        blockStarts = BED.blockstarts(t)
        s_g = BED.chromstart(t) + 1
        L = sum(blockSizes)

        coverage = zeros(Float64, 100)

        # Collect genomic positions on exons
        exon_coordinates = zeros(Int, L)
        offset = 0
        for i in 1:length(blockSizes) # For each exon
            s_g = blockStarts[i] + s_g
            for j in 0:(blockSizes[i]-1)
                exon_coordinates[j+offset+1] = s_g + j
            end
            offset += blockSizes[i]
        end
    end
end

@time unittest3(transcripts[1:100], reader, transcript_length_cut);
# julia> @time unittest3(transcripts[1:100], reader, transcript_length_cut);
# 100
#   0.002450 seconds (1.16 k allocations: 1.887 MiB)

# julia> @time unittest3(transcripts[1:1000], reader, transcript_length_cut);
# 1000
#   0.045619 seconds (11.43 k allocations: 18.670 MiB, 63.43% gc time)


#############
function unittest4(transcripts, reader, transcript_length_cut)
    println(length(transcripts))

    N_bin = 100
	# Define relative gene body coverage output
	relative_coverage = zeros(Float64, N_bin)
    for t in transcripts
        # 
        if sum(BED.blocksizes(t)) < transcript_length_cut
            continue
        end

        if BED.strand(t) != STRAND_POS && BED.strand(t) != STRAND_NEG  # STRAND_NA, STRAND_BOTH
            continue
        end

        blockSizes = BED.blocksizes(t)
        blockStarts = BED.blockstarts(t)
        s_g = BED.chromstart(t) + 1
        L = sum(blockSizes)

        coverage = zeros(Float64, 100)

        # Collect genomic positions on exons
        exon_coordinates = zeros(Int, L)
        offset = 0
        for i in 1:length(blockSizes) # For each exon
            s_g = blockStarts[i] + s_g
            for j in 0:(blockSizes[i]-1)
                exon_coordinates[j+offset+1] = s_g + j
            end
            offset += blockSizes[i]
        end

        # Select percentile positions (resulting in 100 positions)
        percentile_positions = percentile_list(exon_coordinates)

        # Note: Read coverage is saved in `cov` in a 'Left justified' manner,
        #       and left is 5'-end.
        # If the strand of a transcript is '-' STRAND_NEG, read coverage was reversed.
        if BED.strand(t) == STRAND_NEG
            percentile_positions = reverse(percentile_positions) .- 1
        end
        println(length(percentile_positions))
    end
end

@time unittest4(transcripts[1:100], reader, transcript_length_cut);



###############
function unittest5(reader, t)
    coverage_transcript_percentile(reader, t);
end

@time for i in 1:100
     unittest5(reader, transcripts[i])
end
#  14.991470 seconds (538.54 k allocations: 49.682 MiB, 0.07% gc time)


function unittest6(reader, t)
    for i in 1:100
        pos = BED.chromstart(t) + 1
        bamToCoverage_cigarAware_base(reader, BED.chrom(t), pos, pos)[1]
    end
end
@time for i in 1:100
    unittest6(reader, transcripts[i]);
end
#  13.789509 seconds (498.36 k allocations: 43.372 MiB, 0.09% gc time)



function unittest7(reader, t)
    for i in 1:100
        pos = BED.chromstart(t) + 1
        bamToCoverage_cigarAware_position(reader, BED.chrom(t), pos)
    end
end
@time for i in 1:100
    unittest7(reader, transcripts[i]);
end
#  13.938402 seconds (463.90 k allocations: 41.151 MiB, 0.12% gc time)



# bamToCoverage_cigarAware_base() の幅と繰り返し回数はどちらが計算時間にクリティカルか？
function unittest8(reader, t)
    for i in 1:100
        pos = BED.chromstart(t) + 1
        bamToCoverage_cigarAware_base(reader, BED.chrom(t), pos, pos+500)[1]
    end
end
@time for i in 1:100
    unittest8(reader, transcripts[i]);
end
# → あまり変わらない (BAMを読み込みに行くこと自体がクリティカル？)
# @time for i in 1:100
#     unittest8(reader, transcripts[i]);
#     end
# 13.579250 seconds (742.70 k allocations: 103.630 MiB, 0.16% gc time)



function unittest9(reader, t)
    coverage_transcript_percentile2(reader, t);
end

@time for i in 1:100
    println(i)
    unittest9(reader, transcripts[i])
end


tras = load_transcript("/data2/myoshimura/21_Nextflow_test/test2_githubcode_releaseTaskAfter_unstrandedSE/output_190527_nextflow_ramdaQC_SE_unstranded/170202_NS500723_0071_AH33LNBGX2/05_rseqc/gene_bodycoverage/gencode.vM18.primary_assembly.annotation.bed")

function unittest10(tras)
    c = 0
    for t in tras
        if sum(BED.blocksizes(t)) < 100
            continue
        end
        c += 1
    end
    return(c)
end
println(unittest10(tras))
# 134184: julia

# 134100 transcripts finished: python

#############################################
# Coverage sum for a transcript
function unittest11(reader, t)
    sum(coverage_transcript_percentile2(reader, t))
end

covsum = zeros(length(tras));
@time for i in 1:length(tras)
    covsum[i] = unittest11(reader, tras[i])
end

#############################################
# Coverage max for a transcript
function unittest12(reader, t)
    maximum(coverage_transcript_percentile2(reader, t))
end
covmax = zeros(length(tras));
@time for i in 1:length(tras)
    covmax[i] = unittest12(reader, tras[i])
end
# julia> findmax(covmax)
# (28191.0, 122558)

# julia> tras[122558]
# GenomicFeatures.BED.Record:
#    chromosome: chr17
#         start: 39846353
#           end: 39848201
#          name: ENSMUST00000198477.1
#         score: 0
#        strand: -
#   thick start: 39848202
#     thick end: 39848201
#      item RGB: RGB{N0f8}(0.0,0.0,0.0)
#   block count: 1
#       [1]: size=1849, start=1
#
# 
# `/data2/myoshimura/21_Nextflow_test/test2_githubcode_releaseTaskAfter_unstrandedSE/output_190527_nextflow_ramdaQC_SE_unstranded/170202_NS500723_0071_AH33LNBGX2/05_rseqc/gene_bodycoverage/gencode.vM18.primary_assembly.annotation.bed`
# の 122558 行目はrRNA gene である: 
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A39846353%2D39848201&hgsid=792122051_DO9iR6Ls4CmekUGKzeO12tKRAqWA
# RepeatMasker (SSU-rRNA_Hsa )
# 

