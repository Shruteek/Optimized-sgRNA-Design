import csv
from os.path import exists
from GenomeTools import *
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence


def __queryForMetagenomeFile():
    """Repeatedly queries the user for a string filepath, validates each filepath input, takes the first valid filepath,
     and returns its corresponding MetaGenome object."""
    filePathInput = ""
    timesQueried = 0
    while not (exists(filePathInput) and isValidFasta(filePathInput)) and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid metagenome file entered, or file does not exist. Please try again.")
        filePathInput = input("Enter the metagenome file path (e.g. metagen.fasta or /c/Users/metagen.fasta): ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not (exists(filePathInput) and isValidFasta(filePathInput)):
            print("Queried too many times. Exiting...")
    return filePathInput

def __queryForDelimitedGuideFile():
    """Repeatedly queries the user for a string filepath, validates each filepath input, takes the first valid filepath,
    and returns the corresponding filepath."""
    filePathInput = ""
    timesQueried = 0
    while not (exists(filePathInput)) and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid delimiter file entered, or file does not exist. Please try again.")
        filePathInput = input("Enter the delimeter file path (e.g. guidesequences.tsv or /c/Users/sequences.csv): ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not exists(filePathInput):
            print("Queried too many times. Exiting...")
    if not exists(filePathInput):
        return ""
    else:
        return filePathInput


def __queryForGuideSequence(metaGenome):
    """Repeatedly queries the user for a string guide sequence, validates each sequence input, takes the first valid
    sequence, and returns it."""
    guideSequence = ""
    timesQueried = 0
    while not isValidGuideSequence(guideSequence) and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid guide sequence or spacer sequence entered. Please try again.")
        guideSequence = input("Enter a 35 bp guide sequence or a 20 bp spacer sequence: ")
        if len(guideSequence) == 20:
            guideSequence = completeGuideSequence(guideSequence, metaGenome)
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not isValidGuideSequence(guideSequence):
            print("Queried too many times. Exiting...")
    if not isValidGuideSequence(guideSequence):
        return ""
    else:
        return guideSequence


def __queryForOffTargetMethod():
    """Repeatedly queries the user about whether or not to use efficient findOffTarget, validates each boolean input,
     and returns the first valid input."""
    inputChoice = ""
    timesQueried = 0
    while not (inputChoice == "Y" or inputChoice == "N") and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid response given. Please try again.")
        inputChoice = input("Should this program use efficient off-Target searching? Enter Y/N: ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not (inputChoice == "Y" or inputChoice == "N"):
            print("Queried too many times. Exiting...")
    return inputChoice == "Y"


def __queryForDelimiterType():
    """Repeatedly queries the user about which file type the intended spacers file will be, validates each type input,
    and returns the first valid file type."""
    inputType = ""
    timesQueried = 0
    while not (inputType == "CSV" or inputType == "TSV") and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid file type given. Please try again.")
        inputType = input("Please enter a valid delimiter-separated file type (e.g. TSV or CSV): ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not (inputType == "CSV" or inputType == "TSV"):
            print("Queried too many times. Exiting...")
    if inputType == "CSV":
        return ","
    elif inputType == "TSV":
        return "\t"
    else:
        return ","


def __queryForProgramType():
    """Repeatedly queries the user about which program mode they want to execute, if any.
    1: Enter a string guide/spacer sequence plus a metagenome and return its on/off-target scores and heuristics.
    2: Enter a delimiter-separated value format (e.g. CSV or TSV) of guide/spacer sequences plus a metagenome and return
    a table of their on/off-target scores and heuristics.
    3: Exit."""
    exitIndex = 3
    inputMode = "initial"
    timesQueried = 0
    while not (inputMode.isnumeric() and (0 < int(inputMode) <= exitIndex)) and timesQueried < 5:
        print("Below are a list of modes for this program, by their index:")
        print("(1) Enter a string guide/spacer sequence plus a metagenome and return its on/off-target scores and "
              + "heuristics.")
        print("(2) Enter a delimiter-separated value format (e.g. CSV or TSV) of guide/spacer sequences plus a "
              + "metagenome and return a table of their on/off-target scores and heuristics.")
        print("(3) Exit this program.")
        if not timesQueried == 0:
            print("Invalid program type chosen. Please try again.")
        inputMode = input("Please choose which program mode to execute: ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not (inputMode.isnumeric() and (0 < int(inputMode) <= exitIndex)):
            print("Queried too many times. Exiting...")
    if inputMode.isnumeric() and (0 < int(inputMode) <= exitIndex):
        return int(inputMode)


def __queryForSave():
    """Repeatedly queries the user about whether or not to save program data to a file, validates each boolean input,
     and returns the first valid input."""
    inputChoice = ""
    timesQueried = 0
    while not (inputChoice == "Y" or inputChoice == "N") and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid response given. Please try again.")
        inputChoice = input("Should this program save its data to a save file? Enter Y/N: ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not (inputChoice == "Y" or inputChoice == "N"):
            print("Queried too many times. Exiting...")
    return inputChoice == "Y"


def __queryForSaveFile():
    """Repeatedly queries the user about what filepath to save to, validates each filepath input,
    and returns the first valid input."""
    filePathInput = ""
    timesQueried = 0
    validInput = False
    while not validInput and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid save file path entered, or file already exists. Please try again.")
        filePathInput = input("Enter the desired save file path (e.g. saveFile or /c/Users/saveFile), NOT including "
                              + "the type extension: ")
        validInput = not (len(filePathInput) == 0 or exists(filePathInput))
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not validInput:
            print("Queried too many times. Exiting...")
    if not validInput:
        return ""
    else:
        return filePathInput


def __runProgram():
    """Runs the program, from deciding the program mode to returning an output."""
    inputMode = __queryForProgramType()
    if inputMode == 1 or inputMode == 2:
        metaGenPath = __queryForMetagenomeFile()
        metaGen = MetaGenome(metaGenPath, metaGenPath)
        if metaGen.size() == 0:
            print("Empty or invalid metagenome provided. Exiting...")
            return
    elif inputMode == 3:
        print("Exiting program...")
        return
    data = None
    if inputMode == 1:
        guideSequence = __queryForGuideSequence(metaGen)
        if not len(guideSequence) == 35:
            print("No valid guide sequences (or contained spacer sequences) given. Exiting...")
            return
        spacerSequence = SpacerSequence(guideSequence, metaGen)
        data = [[spacerSequence.getGuideSequence(), spacerSequence.getOnTargetScore(),
                 spacerSequence.getOffTargetScores(), spacerSequence.getHeuristic()]]
    elif inputMode == 2:
        delimiter = __queryForDelimiterType()
        delimiterFile = __queryForDelimitedGuideFile()
        if len(delimiterFile) == 0:
            print("No valid delimeter files given. Exiting...")
            return
        read_guides = csv.reader(delimiterFile, delimiter=delimiter)
        guides = []
        for row in read_guides:
            for entry in row:
                guides.append(entry)
    saveToFile = __queryForSave()
    if saveToFile:
        saveFile = __queryForSaveFile()
        if not len(saveFile) == 0:
            print("Program data saved to " + saveToFile(data, saveFile))


if __name__ == '__main__':
    """Handles SeedSequence.py, accepting inputs and outputs to the class and returning the result of on-target and 
    off-target calculations and heuristics."""
    __runProgram()
