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
    elif not (len(coreSequence) == 35 and coreSequence[27:29] == "CC"):
        print("Given RNA sequence is not 35 bp with a PAM: " + coreSequence)
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


def crossCorrelateSequences(spacerSequence, substrateSequence):
    """Takes in a DNA string spacerSequence and a DNA string substrateSequence, cross-correlates keySequence along
        substrateSequence, and returns the numerical results as a list of integers."""
    if not isValidDNA(spacerSequence) or not isValidDNA(substrateSequence):
        print("Given sequences are not valid DNA: " + spacerSequence + " and " + substrateSequence)
        return []
    elif len(spacerSequence) != 20:
        print("Given spacerSequence is not 20 nucleotides long.")
    else:
        crosscorrelation = []
        for shift in range(0, len(substrateSequence) - 20 - 3 - 6 - 6 + 1):
            if substrateSequence[(shift + 27):(shift + 29)] == "GG":
                correlation = correlateSequences(spacerSequence, substrateSequence[(shift + 6):(shift + 26)])
                crosscorrelation.append(correlation)
        return crosscorrelation


def correlateSequences(keySequence, substrateSubsequence):
    """Takes in an arbitrary string keySequence and an arbitrary string substrate and calculates the correlation
    between them, where equivalent characters add 1 and other characters add 0, then returns the correlation."""
    correlation = 0
    for nucleotideIndex in range(len(keySequence)):
        if keySequence[nucleotideIndex] == substrateSubsequence[nucleotideIndex]:
            correlation = correlation + 1
        elif substrateSubsequence[nucleotideIndex] == "N" or keySequence[nucleotideIndex] == "N":
            correlation = correlation + 1
    return correlation
