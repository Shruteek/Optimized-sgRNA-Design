import csv
import re
import os


def isValidFasta(FASTAFile):
    """Takes in a string filename and returns whether the given path points to an existing FASTA file."""
    return os.path.exists(FASTAFile) \
           and str.lower(os.path.splitext(FASTAFile)[1]) == ".fasta" \
           or str.lower(os.path.splitext(FASTAFile)[1]) == ".fastq" \
           or str.lower(os.path.splitext(FASTAFile)[1]) == ".gz"


def isValidCSV(CSVFile):
    """Takes in a string filename and returns whether the given path points to an existing CSV."""
    return os.path.exists(CSVFile) and str.lower(os.path.splitext(CSVFile)[1]) == ".csv"


def isValidTSV(TSVFile):
    """Takes in a string filename and returns whether the given path points to an existing TSV."""
    return os.path.exists(TSVFile) and str.lower(os.path.splitext(TSVFile)[1]) == ".tsv"


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


def isValidTargetSpacerInput(inputSequence):
    """Method that returns whether a given sequence is a valid string representing any of the following:
    1. A target/spacer sequence (20 base-pair RNA/DNA sequence resembling the target/spacer sequence)
    2. A target/spacer + PAM sequence (23 base-pair RNA/DNA sequence resembling the target sequence with its PAM)
    Returns True if the given sequence is any of the above and False otherwise."""
    if isValidDNA(inputSequence) or isValidRNA(inputSequence):
        if len(inputSequence) == 20:
            return True
        elif len(inputSequence) == 23:
            return inputSequence[21:23] == "GG"
    return False


def complementaryDNA(sequence):
    """Method that returns the complementary DNA strand sequence to the given input DNA or RNA sequence."""
    complement = ""
    if not (isValidRNA(sequence) or isValidDNA(sequence)):
        return complement
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


def reverseComplementaryDNA(sequence):
    """Method that returns the reverse of the complementary DNA strand sequence to the given input DNA or RNA
    sequence. """
    reverseComplement = ""
    if not (isValidRNA(sequence) or isValidDNA(sequence)):
        return reverseComplement
    for nucleotide in str.upper(sequence):
        if nucleotide == "A":
            reverseComplement = reverseComplement + "T"
        elif nucleotide == "C":
            reverseComplement = reverseComplement + "G"
        elif nucleotide == "T":
            reverseComplement = reverseComplement + "A"
        elif nucleotide == "U":
            reverseComplement = reverseComplement + "A"
        elif nucleotide == "G":
            reverseComplement = reverseComplement + "C"
    return reverseComplement[::-1]


def convertToDNA(sequence):
    """Method that converts the given input RNA or DNA sequence to DNA."""
    complement = ""
    if not (isValidRNA(sequence) or isValidDNA(sequence)):
        return complement
    for nucleotide in str.upper(sequence):
        if nucleotide == "A":
            complement = complement + "A"
        elif nucleotide == "C":
            complement = complement + "C"
        elif nucleotide == "U":
            complement = complement + "T"
        elif nucleotide == "T":
            complement = complement + "T"
        elif nucleotide == "G":
            complement = complement + "G"
    return complement


def appendSpacerToData(data, spacer):
    """Takes in a data list and a SpacerSequence object, and appends the relevant data from the object to the list
    in a tidy data format, where each row represents an on/off-target sequence and its relevant information."""
    if not isinstance(data, list):
        print("Incorrect type given for data: " + type(data))
        return data
    data.append(["Spacer",
                 "Total Alignments",
                 "Sequence_Type",
                 "Location",
                 "Sequence",
                 "Surrounding_Sequence",
                 "Mismatches",
                 "On_Target_Score",
                 "Off_Target_Score"])
    targets = spacer.getTargets()
    spacer_sequence = spacer.getSpacerSequence()
    total_alignments = spacer.getTotalAlignments()
    for targetSequence in targets:
        data.append([spacer_sequence,
                     total_alignments,
                     targetSequence.sequence_type,
                     targetSequence.location,
                     targetSequence.target_sequence,
                     targetSequence.surrounding_sequence,
                     targetSequence.mismatches,
                     targetSequence.on_target_score,
                     targetSequence.off_target_score])
    return data


def writeNestedListToRows(nestedList, saveFilePath):
    """Takes in a nested list containing data and a save file name (WITH the file type extension), writes the
    nested list data to a file at the given saveFilePath (after the contents if the file already exists),
    and returns the path of the saved file."""
    if isValidCSV(saveFilePath):
        writeNestedListToCSVRows(nestedList, saveFilePath)
    elif isValidTSV(saveFilePath):
        writeNestedListToTSVRows(nestedList, saveFilePath)
    else:
        print("Invalid TSV/CSV file input.")


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

