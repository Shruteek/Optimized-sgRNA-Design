import csv
import re
from os.path import exists


def isValidFasta(FASTAFile):
    """Takes in a string filename and returns whether the given path points to an existing FASTA file."""
    return isinstance(FASTAFile, str) and exists(FASTAFile)


def isValidCSV(CSVFile):
    """Takes in a string filename and returns whether the given path points to an existing CSV."""
    return isinstance(CSVFile, str) and len(CSVFile) > 4 and str.lower(CSVFile[-4:]) == ".csv" and exists(CSVFile)


def isValidTSV(TSVFile):
    """Takes in a string filename and returns whether the given path points to an existing TSV."""
    return isinstance(TSVFile, str) and len(TSVFile) > 4 and str.lower(TSVFile[-4:]) == ".tsv" and exists(TSVFile)


def isValidDNA(DNASequence):
    """Method that returns a boolean representing whether the input sequence is a valid DNA sequence."""
    if not isinstance(DNASequence, str):
        return False
    validDNABasePairs = "ACTGN"
    for nucleotideIndex in range(len(DNASequence)):
        if not validDNABasePairs.__contains__(DNASequence[nucleotideIndex]):
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


def isValidGuideSequence(guideSequence):
    """Method that returns whether a given sequence is a valid string representing a 35 base-pair RNA sequence in the
        format 6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream"""
    if not isinstance(guideSequence, str):
        # print("Given sequence is not a valid string: " + coreSequence)
        return False
    elif not isValidRNA(guideSequence):
        # print("Given guide sequence is not valid RNA: " + coreSequence)
        return False
    elif not (len(guideSequence) == 35 and guideSequence[27:29] == "CC"):
        # print("Given RNA sequence is not 35 bp with a PAM: " + coreSequence)
        return False
    else:
        return True


def complementaryRNA(sequence):
    """Method that returns the complementary RNA strand sequence to the given input RNA or DNA sequence."""
    complement = ""
    if not (isValidRNA(sequence) or isValidDNA(sequence)):
        return ""
    for nucleotide in str.upper(sequence):
        if nucleotide == "A":
            complement = complement + "U"
        elif nucleotide == "C":
            complement = complement + "G"
        elif nucleotide == "U":
            complement = complement + "A"
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
    """Takes in a 20-bp RNA spacer sequence and a metaGenome and returns the full RNA guide sequence according to the
    metaGenome, if the complement of the spacer exists in the metagenome and is within usable bounds."""
    if not isValidRNA(spacer):
        print("Spacer is not a string.")
        return spacer
    elif not len(spacer) == 20:
        print("Spacer is not 20 base pairs.")
        return spacer
    else:
        target = complementaryDNA(spacer)
        matchingSequences = metaGenome.getSequence(target)
        matchingSubsequenceIndices = []
        for sequence in matchingSequences:
            matchingSubsequenceIndices.append([m.start() for m in re.finditer(target, sequence)])
        for sequenceIndex in range(len(matchingSequences)):
            for subsequenceIndex in matchingSubsequenceIndices[sequenceIndex]:
                if 6 < subsequenceIndex < len(matchingSequences[sequenceIndex]) - 6 and \
                        matchingSequences[sequenceIndex][subsequenceIndex:(subsequenceIndex+20)] == target and \
                        matchingSequences[sequenceIndex][(subsequenceIndex+21):(subsequenceIndex+20+3)] == "GG":
                    return complementaryRNA(matchingSequences[sequenceIndex][(subsequenceIndex-6):(subsequenceIndex+29)])
        return spacer


def findTargetsFromSpacer(spacerSequence, substrateSequence, mismatchStrictness):
    """Takes in a 20 bp RNA string spacerSequence and a substrateSequence, cross-correlates the DNA complement of the
    spacer with the full substrate, finds each 35 bp subsequence that has a PAM (NGG) and a  correlation of at least
    20 - mismatchStrictness, and returns a list of valid 35 bp target guide sequences."""
    if not (isValidRNA(spacerSequence) and len(spacerSequence) == 20):
        # print("Given sequence is not a valid spacer sequence: " + spacerSequence)
        return []
    else:
        targetSequence = complementaryDNA(spacerSequence)
        offTargets = []
        crosscorrelation = crossCorrelateSequencesEfficiently(targetSequence, substrateSequence)
        for shift in range(len(crosscorrelation)):
            if crosscorrelation[shift] >= len(targetSequence) - mismatchStrictness:
                # print("Off-target sequence: " + substrateSequence[shift + 6:(shift + 6 + 20 + 3)]
                #       + " against " + targetSequence)
                # print("with a correlation of " + str(crosscorrelation[shift]))
                offTargets.append(substrateSequence[shift:(shift + 6 + 6 + 20 + 3)])
        return offTargets


def crossCorrelateSequences(targetSequence, substrateSequence):
    """Takes in a 20 bp DNA string targetSequence and an arbitrary length DNA string substrateSequence, cross-correlates
    targetSequence by shifting it along substrateSequence, and returns the results. Correlates the 20 bp ideal target
    sequence with the 20 bp in the middle of every potential 35 bp target guide sequence on substrateSequence, and
    returns the numerical correlation results as a list of integers, with each correlation result between 20 bp
    targetSequence and substrate 35 bp subsequence being stored at the cross-correlation shift index of the list.
    Records 0 for any 35 bp region without a PAM."""
    if not isValidDNA(targetSequence) or not isValidDNA(substrateSequence):
        # print("Given sequences are not valid DNA: " + targetSequence + " and " + substrateSequence)
        return []
    elif len(targetSequence) != 20:
        # print("Given targetSequence is not 20 nucleotides long.")
        return []
    else:
        crosscorrelation = []
        for shift in range(0, len(substrateSequence) - 20 - 3 - 6 - 6 + 1):
            correlation = 0
            if substrateSequence[(shift + 27):(shift + 29)] == "GG":
                correlation = correlateSequences(targetSequence, substrateSequence[(shift + 6):(shift + 26)])
            crosscorrelation.append(correlation)
        return crosscorrelation


def crossCorrelateSequencesEfficiently(targetSequence, substrateSequence):
    """Takes in a 20 bp DNA string targetSequence and an arbitrary length DNA string substrateSequence,
    cross-correlates targetSequence by shifting it along substrateSequence, and returns the results. Correlates the
    first PAMAdjacentLength bp in the 20 bp ideal target sequence (AKA the PAM-adjacent region) with the first
    PAMAdjacentLength bp of the 20 bp in the middle of every potential 35 bp target guide  sequence on
    substrateSequence, and if there is no more than PAMMismatchStrictness mismatches, returns the full numerical
    correlation of both 20 bp regions (target sequence and potential target guide subsequence), storing the results
    as a list of integers, with each correlation result between 20 bp targetSequence and substrate 35 bp subsequence
    being stored at the cross-correlation shift index of the list. Records 0 for any 35 bp region without a PAM.
    Records just the 10 bp PAM-adjacent correlation for any region with more than 1 PAM-adjacent mismatch. """
    if not isValidDNA(targetSequence):
        return []
    elif not isValidDNA(substrateSequence):
        return []
    elif len(targetSequence) != 20:
        return []
    else:
        crosscorrelation = []
        PAMAdjacentLength = 10
        PAMMismatchStrictness = 1
        for shift in range(0, len(substrateSequence) - 20 - 3 - 6 - 6 + 1):
            correlation = 0
            if substrateSequence[(shift + 27):(shift + 29)] == "GG":
                correlation = correlateSequences(targetSequence[0:PAMAdjacentLength],
                                                 substrateSequence[(shift + 6):(shift + 6 + PAMAdjacentLength)])
                if correlation >= PAMAdjacentLength - PAMMismatchStrictness:
                    correlation = correlateSequences(targetSequence, substrateSequence[(shift + 6):(shift + 26)])
            crosscorrelation.append(correlation)
        return crosscorrelation


def correlateSequences(targetSequence, substrateSubsequence):
    """Takes in any two DNA strings (targetSequence and substrateSubsequence) of the same length, validates the inputs,
    and calculates the correlation between them, where equivalent characters add 1 and other characters add 0, then
    returns the correlation of the two strings."""
    correlation = 0
    if not (isinstance(targetSequence, str) and isinstance(substrateSubsequence, str)):
        # print("Given sequences to correlate are not strings.")
        return correlation
    elif not len(targetSequence) == len(substrateSubsequence):
        # print("Sequences to correlate are not same length: " + targetSequence + " and " + substrateSubsequence)
        return correlation
    else:
        for nucleotideIndex in range(len(targetSequence)):
            if targetSequence[nucleotideIndex] == substrateSubsequence[nucleotideIndex]:
                correlation = correlation + 1
            elif substrateSubsequence[nucleotideIndex] == "N" or targetSequence[nucleotideIndex] == "N":
                correlation = correlation + 1
    return correlation


def saveToTXT(data, saveFilePath):
    """Takes in arbitrary data and a save file name (WITHOUT the file type extension), writes the data to a .txt file at
    the given saveFilePath, and returns the path of the saved .txt file."""
    saveFile = open(saveFilePath + ".txt", "w+")
    saveFile.write(data)
    return saveFilePath + ".txt"


def writeNestedListToCSVRows(nestedList, saveFilePath):
    """Takes in a nested list containing data and a .CSV save file name (WITH the file type extension), writes the
    nested list data to a .CSV file at the given saveFilePath (after the contents if the file already exists),
    and returns the path of the saved .CSV file."""
    CSVFile = open(saveFilePath, "a", newline="")
    CSVlinewriter = csv.writer(CSVFile, delimiter=",")
    for sublist in nestedList:
        CSVlinewriter.writerow(sublist)
    CSVFile.close()
    return saveFilePath


def writeNestedListToTSVRows(nestedList, saveFilePath):
    """Takes in a nested list containing data and a .TSV save file name (WITH the file type extension), writes the
    nested list data to a .TSV file at the given saveFilePath (after the contents if the file already exists),
    and returns the path of the saved .TSV file."""
    TSVFile = open(saveFilePath, "a+", newline="")
    TSVlinewriter = csv.writer(TSVFile, delimiter=" ")
    for sublist in nestedList:
        TSVlinewriter.writerow(sublist)
    TSVFile.close()
    return saveFilePath


def writeListToCSVRow(contentsList, saveFilePath):
    """Takes in a list containing data and a .CSV save file name (WITH the file type extension), appends each element
        of the list to a cell in a new row of a .CSV file at the given saveFilePath (keeping the contents if the file
        already exists), and returns the path of the saved .CSV file."""
    TSVFile = open(saveFilePath, "a+", newline="")
    TSVlinewriter = csv.writer(TSVFile, delimiter=",")
    TSVlinewriter.writerow(contentsList)
    TSVFile.close()
    return saveFilePath


def writeListToTSVRow(contentsList, saveFilePath):
    """Takes in a list containing data and a .TSV save file name (WITH the file type extension), appends each element
        of the list to a cell in a new row of a .TSV file at the given saveFilePath (keeping the contents if the file
        already exists), and returns the path of the saved .TSV file."""
    TSVFile = open(saveFilePath, "a", newline="")
    TSVlinewriter = csv.writer(TSVFile, delimiter=" ")
    TSVlinewriter.writerow(contentsList)
    TSVFile.close()
    return saveFilePath
