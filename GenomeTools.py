def isValidFasta(sequencesFile):
    """Takes in a string filename and returns whether the given variable is a valid fasta filename."""
    return isinstance(sequencesFile, str) and ".fasta" in str.lower(sequencesFile) and \
           str.lower(sequencesFile[-6:]) == ".fasta" and len(sequencesFile) > 6


def isValidDNA(DNASequence):
    """Method that returns a boolean representing whether the input sequence is a valid DNA sequence."""
    if not isinstance(DNASequence, str):
        return False
    validDNABasePairs = "ACTGN"
    for nucleotide in DNASequence:
        if not validDNABasePairs.__contains__(nucleotide):
            return False
    return True


def isValidRNA(RNASequence):
    """Method that returns a boolean representing whether the input sequence is a valid RNA sequence."""
    if not isinstance(RNASequence, str):
        return False
    validRNABasePairs = "ACUG"
    for nucleotide in RNASequence:
        if not validRNABasePairs.__contains__(nucleotide):
            return False
    return True


def isValidCoreSequence(coreSequence):
    """Method that returns whether a given sequence is a valid string representing a 35 base-pair RNA sequence in the
        format 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream"""
    if not isinstance(coreSequence, str):
        print("Given sequence is not a valid string: " + coreSequence)
        return False
    elif not isValidRNA(coreSequence):
        print("Given sequence is not valid RNA: " + coreSequence)
    elif not (len(coreSequence) == 35 and coreSequence[27:28] == "GG"):
        print("Given RNA sequence is not valid RNA: " + coreSequence)
        return False
    else:
        return True


def complementaryRNA(sequence):
    """Method that returns the complementary RNA strand sequence to the given input DNA sequence."""
    complement = ""
    for nucleotide in str.upper(sequence):
        if nucleotide == "A":
            complement = complement + "U"
        elif nucleotide == "C":
            complement = complement + "G"
        elif nucleotide == "T":
            complement = complement + "A"
        elif nucleotide == "G":
            complement = complement + "C"
    return complement


def complementaryDNA(sequence):
    """Method that returns the complementary DNA strand sequence to the given input DNA or RNA sequence."""
    complement = ""
    for nucleotide in str.upper(sequence):
        if nucleotide == "A":
            complement = complement + "T"
        elif nucleotide == "C":
            complement = complement + "G"
        elif nucleotide == "T":
            complement = complement + "A"
        elif nucleotide == "U":
            complement = complement + "A"
        elif nucleotide == "G":
            complement = complement + "C"
    return complement


def correlateSequences(keySequence, substrateSequence):
    """Takes in a string keySequence and a string substrateSequence, cross-correlates keySequence along
        substrateSequence, and returns the numerical results as a list of integers."""
    if not isValidDNA(None, keySequence) or not isValidDNA(None, substrateSequence):
        print("Given sequences are not valid DNA: " + keySequence + " and " + substrateSequence)
        return []
    else:
        crosscorrelation = []
        for shift in range(len(substrateSequence) - len(keySequence) + 1):
            correlation = 0
            for nucleotideIndex in range(len(keySequence)):
                if keySequence[nucleotideIndex] == substrateSequence[nucleotideIndex + shift]:
                    correlation = correlation + 1
                elif substrateSequence[nucleotideIndex + shift] == "n":
                    correlation = correlation + 1
            crosscorrelation.append(correlation)
        return crosscorrelation