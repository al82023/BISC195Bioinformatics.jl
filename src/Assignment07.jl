module Assignment07

export normalizeDNA,
        composition,
        gc_content,
        complement,
        reverse_complement

# # uncomment the following line if you intend to use BioSequences types
using BioSequences
import BioSequences: composition, gc_content, complement, reverse_complement

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
    BioSequences.composition(normalizeDNA(sequence))
end

function gc_content(sequence)
    BioSequences.gc_content(normalizeDNA(sequence))
end

function complement(sequence::AbstractString)
    BioSequences.complement(normalizeDNA(sequence))
end

function reverse_complement(sequence::AbstractString)
    BioSequences.reverse_complement(normalizeDNA(sequence))
end

end # module Assignment07