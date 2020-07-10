module Assignment07

export normalizeDNA
        composition

# # uncomment the following line if you intend to use BioSequences types
using BioSequences
import BioSequences: composition

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


function composition(sequence)
    BioSequences.composition(normalizeDNA(sequence))
end

end # module Assignment07