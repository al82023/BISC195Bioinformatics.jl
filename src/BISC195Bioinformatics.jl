module BISC195Bioinformatics

export normalizeDNA,
        composition,
        gc_content,
        complement,
        reverse_complement,
        parse_fasta,
        LongDNASeq,
        uniquekmers,
        distance

using BioSequences
import BioSequences: composition, gc_content, complement, reverse_complement, LongDNASeq

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return LongDNASeq(seq) # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end


# Your code here.
# Don't forget to export your functions!


function composition(sequence::AbstractString)
    BioSequences.composition(LongDNASeq(sequence))
end

function gc_content(sequence::AbstractString)
    BioSequences.gc_content(LongDNASeq(sequence))
end

function complement(sequence::AbstractString)
    BioSequences.complement(LongDNASeq(sequence))
end

function complement(base::Char)
    nt = convert(DNA, base)
    Char(BioSequences.complement(nt))
end

function reverse_complement(sequence::AbstractString)
    BioSequences.reverse_complement(LongDNASeq(sequence))
end

function parse_fasta(path)
    headers = []
    sequences = []
    str = ""
    for line in eachline(path)
        if '>' ∈ line
            push!(headers, string(lstrip(line, '>')))
            if str != ""
                push!(sequences, LongDNASeq(str))
            end
            str = ""
        elseif line != ""
            # for base in line
            #     occursin(base, "AGCTNagctn") || error("invalid base $base")
            # end
            str = str * line
        end
    end
    push!(sequences, LongDNASeq(str))
    parsedfile = (headers, sequences)
    return parsedfile
end

function uniquekmers(seq, k)
    sequence = String(seq)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = []
    stopindex = length(sequence)-(k-1)
    for i in 1:stopindex
        kmer = sequence[i:i+k-1]
        if match(r"[^ACGT]", kmer) != nothing # if kmer contains any ambiguous base
            continue
        end
        push!(kmers, LongDNASeq(kmer))
    end
    return Set(kmers)
end

function distance(set1, set2)
    1 - (length(intersect(set1, set2)) / length(union(set1, set2)))
end

end # module BISC195Bioinformaticsp