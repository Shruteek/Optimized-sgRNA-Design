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


def isValidGuideSequence(coreSequence):
    """Method that returns whether a given sequence is a valid string representing a 35 base-pair RNA sequence in the
        format 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream"""
    if not isinstance(coreSequence, str):
        # print("Given sequence is not a valid string: " + coreSequence)
        return False
    elif not isValidRNA(coreSequence):
        # print("Given guide sequence is not valid RNA: " + coreSequence)
        return False
    elif not (len(coreSequence) == 35 and coreSequence[27:29] == "CC"):
        # print("Given RNA sequence is not 35 bp with a PAM: " + coreSequence)
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


def completeGuideSequence(spacer, metaGenome):
    """Takes in a 20-bp spacer sequence and a metaGenome and returns the full RNA guide sequence according to the
    metaGenome, if the complement of the spacer exists in the metagenome and is within usable bounds."""
    if not isinstance(spacer, str):
        # print("Spacer is not a string.")
        return spacer
    elif not len(spacer) == 20:
        # print("Spacer is not 20 base pairs.")
        return spacer
    else:
        matchingSequences = metaGenome.getSequence(complementaryDNA(spacer))
        if len(matchingSequences) == 0:
            return spacer
        else:
            guideIndex = matchingSequences[0].find(complementaryDNA(spacer)) - 6
            if guideIndex < 0 or guideIndex > len(matchingSequences[0]) - 35:
                return spacer
            else:
                return complementaryRNA(matchingSequences[0][guideIndex:(guideIndex + 35)])


def crossCorrelateSequences(targetSequence, substrateSequence):
    """Takes in a DNA string targetSequence and a DNA string substrateSequence, cross-correlates keySequence along
        substrateSequence, and returns the numerical results as a list of integers."""
    if not isValidDNA(targetSequence) or not isValidDNA(substrateSequence):
        print("Given sequences are not valid DNA: " + targetSequence + " and " + substrateSequence)
        return []
    elif len(targetSequence) != 20:
        print("Given targetSequence is not 20 nucleotides long.")
        return []
    else:
        crosscorrelation = []
        for shift in range(0, len(substrateSequence) - 20 - 3 - 6 - 6 + 1):
            correlation = 0
            if substrateSequence[(shift + 27):(shift + 29)] == "GG":
                correlation = correlateSequences(targetSequence, substrateSequence[(shift + 6):(shift + 26)])
            crosscorrelation.append(correlation)
        return crosscorrelation


def correlateSequences(keySequence, substrateSubsequence):
    """Takes in an arbitrary string keySequence and an arbitrary string substrate and calculates the correlation
    between them, where equivalent characters add 1 and other characters add 0, then returns the correlation."""
    correlation = 0
    if not (isinstance(keySequence, str) and isinstance(substrateSubsequence, str)):
        print("Given sequences to correlate are not strings.")
    elif not len(keySequence) == len(substrateSubsequence):
        print("Sequences to correlate are not same length: " + keySequence + " and " + substrateSubsequence)
    else:
        for nucleotideIndex in range(len(keySequence)):
            if keySequence[nucleotideIndex] == substrateSubsequence[nucleotideIndex]:
                correlation = correlation + 1
            elif substrateSubsequence[nucleotideIndex] == "N" or keySequence[nucleotideIndex] == "N":
                correlation = correlation + 1
    return correlation


def saveToFile(data, saveFilePath):
    """Takes in arbitrary data and a save file name (WITHOUT the file type extension), writes the data to a .txt file at
    the given saveFilePath, and returns the path of the saved .txt file."""
    saveFile = open(saveFilePath + ".txt", "w+")
    saveFile.write(data)
    return saveFilePath + ".txt"