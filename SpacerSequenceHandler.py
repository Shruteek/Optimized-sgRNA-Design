import sys
import os
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
        filePathInput = input("Enter the delimiter file path (e.g. guidesequences.tsv or /c/Users/sequences.csv): ")
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
        guideSequence = input("Enter a 20 bp spacer sequence: ")
        if len(guideSequence) == 20:
            guideSequence = completeGuideSequence(guideSequence, metaGenome)
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not isValidGuideSequence(guideSequence):
            print("Queried too many times. Exiting...")
    if not isValidGuideSequence(str.upper(guideSequence)):
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
        filePathInput = input("Enter the desired save file path (e.g. saveFile.CSV or /c/Users/saveFile.TSV), "
                              "including the type extension: ")
        validInput = not (isValidCSV(filePathInput) or isValidTSV(filePathInput))
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not validInput:
            print("Queried too many times. Exiting...")
    if not validInput:
        return ""
    else:
        return filePathInput


def __printHelpMessage(helpArgument):
    """Simply prints out one of the help messages explaining program details, returning nothing."""
    if helpArgument == "no-input":
        print("No input (e.g. `python SpacerSequenceHandler.py`):\n"
              "In no input mode, the program can take a 20 nucleotide character string representing the 20 base pair "
              "sgRNA spacer sequence and a metagenome .FASTA file and return the target analysis list of the sgRNA. It "
              "can also take a delimiter-type file (e.g. CSV or TSV) of sgRNA spacers and a metagenome .FASTA file and "
              "return a list of target analysis lists, one for each valid sgRNA. These inputs are repeatedly queried "
              "for command line text input until 5 invalid inputs are given for a prompt, or until all valid inputs "
              "have been given.")
    elif helpArgument == "single-spacer-input" or helpArgument == "SSI":
        print("Single spacer input (e.g `python SpacerSequenceHandler.PY SSI ACTGACTGACTGACTGACTG metaGenome.FASTA "
              "saveFile.CSV`):\n"
              "In single spacer input mode, the program can take a 20 nucleotide character string representing the "
              "20 base pair sgRNA spacer sequence (though currently, no generic 'N' functionality), a .FASTA file path "
              "representing the metaGenome with which the spacer sequence is to be analyzed, and (optionally) a .CSV "
              "or .TSV file path representing the file to which the target analysis data should be written (**if the "
              "save file already exists, its contents will be overwritten**). The target analysis data will take the "
              "form:\n"
              "[35 bp RNA guide sequence (6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream); "
              "on-target score; list of off-target scores; RNA guide sequence heuristic]")
    elif helpArgument == "multi-spacer-input" or helpArgument == "MSI":
        print("Multi spacer input (e.g `python SpacerSequenceHandler.PY MSI spacerSequences.CSV metaGenome.FASTA "
              "saveFile.CSV`):\n"
              "In multi spacer input mode, the program can take a .TSV or .CSV file path representing"
              "nucleotide characters (though currently, no functionality for the 'N' character), a .FASTA file path "
              "representing the metaGenome with which the spacer sequences are to be analyzed, and (optionally) a .CSV "
              "or .TSV file path representing the file to which the target analysis data should be written (**if the "
              "save file already exists, its contents will be overwritten**). Each target analysis row will take the "
              "form:\n"
              "[35 bp RNA guide sequence (6 bp upstream, 20 bp spacer sequence, 3 bp PAM, 6 bp downstream); "
              "on-target score; list of off-target scores; RNA guide sequence heuristic]")
    else:
        if not helpArgument == "":
            print("[Error] Invalid help argument Entered. See generic help message below.")
        print("This program can be run with any of the below input schemes.\n"
              "No input (e.g. `python SpacerSequenceHandler.PY`): The program will repeatedly prompt for inputs until "
              "it receives all necessary information or too many failed inputs.\n\n"
              "Single spacer input (e.g `python SpacerSequenceHandler.PY SSI ACTGACTGACTGACTGACTG metaGenome.FASTA "
              "saveFile.CSV`): The program takes in a 20 bp sgRNA spacer and a .FASTA file path and saves the target "
              "analysis data to the saveFile path.\n\n "
              "Multi spacer input (e.g `python SpacerSequenceHandler.PY MSI spacerSequences.CSV metaGenome.FASTA "
              "saveFile.CSV`): The program takes in a .CSV or .TSV file path and a .FASTA file path and saves a list "
              "of target analysis data to the saveFile path.\n\n"
              "Help (e.g. `python SpacerSequenceHandler.py help`): Print this message. For more detailed help on "
              "specific modes, use one of the following commands:\n"
              "`python SpacerSequenceHandler.py help no-input`\n"
              "`python SpacerSequenceHandler.py help single-spacer-input`\n"
              "`python SpacerSequenceHandler.py help multi-spacer-input")


def __runProgram():
    """Runs the program, from deciding the program mode to returning an output."""
    inputMode = __queryForProgramType()
    if inputMode == 1 or inputMode == 2:
        metaGenPath = __queryForMetagenomeFile()
        metaGenome = MetaGenome(metaGenPath, metaGenPath)
        data = []
        if metaGenome.size() == 0:
            print("Empty or invalid metagenome provided. Exiting...")
            return
        elif inputMode == 1:
            guideSequence = __queryForGuideSequence(metaGenome)
            if not len(guideSequence) == 35:
                print("No valid guide sequences (or contained spacer sequences) given. Exiting...")
                return
            spacerSequence = SpacerSequence(guideSequence, metaGenome)
            data.append([spacerSequence.getGuideSequence(), spacerSequence.getOnTargetScore(),
                         spacerSequence.getOffTargetScores(), spacerSequence.getHeuristic()])
        elif inputMode == 2:
            delimiter = __queryForDelimiterType()
            delimiterFile = __queryForDelimitedGuideFile()
            if len(delimiterFile) == 0:
                print("No valid delimiter files given. Exiting...")
                return
            guidesTable = csv.reader(open(delimiterFile, newline=''), delimiter=delimiter)
            rownum = 0
            for guideRow in guidesTable:
                rownum = rownum + 1
                columnnum = 0
                for guideColumnEntry in guideRow:
                    columnnum = columnnum + 1
                    if isValidGuideSequence(guideColumnEntry):
                        print("GuideEntry [" + str(rownum) + ", " + str(columnnum) + "]: " + guideColumnEntry)
                        spacerSequenceColumnEntry = SpacerSequence(guideColumnEntry, metaGenome)
                        data.append([spacerSequenceColumnEntry.getGuideSequence(),
                                     spacerSequenceColumnEntry.getOnTargetScore(),
                                     spacerSequenceColumnEntry.getOffTargetScores(),
                                     spacerSequenceColumnEntry.getHeuristic()])
        doSaveToFile = __queryForSave()
        if doSaveToFile:
            saveFilePath = __queryForSaveFile()
            if not len(saveFilePath) == 0:
                if isValidTSV(saveFilePath):
                    saveNestedListToTSV(data, saveFilePath)
                elif isValidCSV(saveFilePath):
                    saveNestedListToCSV(data, saveFilePath)
                print("Program data saved to " + saveFilePath)
    elif inputMode == 3:
        print("Exiting program...")
        return


if __name__ == '__main__':
    """Handles SeedSequence.py, accepting inputs and outputs to the class and returning the result of on-target and 
    off-target calculations and heuristics."""
    # print('Number of arguments:', len(sys.argv), 'arguments.')
    # print('Argument List:', str(sys.argv))
    if len(sys.argv) == 1:
        __runProgram()
    elif len(sys.argv) >= 2:
        if sys.argv[1] == "help":
            if len(sys.argv) == 2:
                __printHelpMessage("")
            elif len(sys.argv) == 3:
                __printHelpMessage(sys.argv[2])
            elif len(sys.argv) > 3:
                print("[Error] Unexpected arguments after help argument: " + str(sys.argv[3:]))
        elif sys.argv[1] == "SSI" or sys.argv[1] == "MSI":
            if len(sys.argv) == 2:
                print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
            elif len(sys.argv) == 3:
                print("[Error] Missing metagenome file input.")
            elif len(sys.argv) > 5:
                print("[Error] Unexpected arguments after metagenome file input: " + str(sys.argv[5:]))
            elif sys.argv[1] == "SSI":
                if not (isValidRNA(sys.argv[2]) and (len(sys.argv[2]) == 20 or len(sys.argv[2]) == 35)):
                    print("[Error] Invalid sgRNA sequence given: " + sys.argv[2] +
                          " (expected 20 or 35 RNA base pairs)")
                elif not isValidFasta(sys.argv[3]):
                    print("[Error] Invalid FASTA file path given: " + sys.argv[3])
                else:
                    metaGen = MetaGenome(sys.argv[3])
                    spacerSequence = SpacerSequence(sys.argv[2], metaGen)
                    data = [[spacerSequence.getGuideSequence(),
                             spacerSequence.getOnTargetScore(),
                             spacerSequence.getOffTargetScores(),
                             spacerSequence.getHeuristic()]]
                    if len(sys.argv) == 5:
                        saveFile = open(sys.argv[4], "w")
                        if isValidCSV(sys.argv[4]):
                            saveFile.close()
                            saveNestedListToCSV(data, sys.argv[4])
                            print(data)
                            print("Successfully saved to " + sys.argv[4])
                        elif isValidTSV(sys.argv[4]):
                            saveFile.close()
                            saveNestedListToTSV(data, sys.argv[4])
                            print(data)
                            print("Successfully saved to " + sys.argv[4])
                        else:
                            saveFile.close()
                            os.remove(sys.argv[4])
                            print(data)
                            print("[Error] Could not save data. Invalid TSV or CSV file path: " + sys.argv[4])
                    else:
                        print(data)
            elif sys.argv[1] == "MSI":
                if not (isValidCSV(sys.argv[2]) or isValidTSV(sys.argv[2])):
                    print("[Error] Invalid sgRNA CSV/TSV file path given: " + sys.argv[2])
                elif not isValidFasta(sys.argv[3]):
                    print("[Error] Invalid FASTA file path given: " + sys.argv[3])
                else:
                    data = []
                    metaGen = MetaGenome(sys.argv[3])
                    if isValidCSV(sys.argv[2]):
                        read_guides = csv.reader(open(sys.argv[2], newline=''), delimiter=',')
                    else:
                        read_guides = csv.reader(open(sys.argv[2], newline=''), delimiter=' ')
                    for row in read_guides:
                        for guideEntry in row:
                            spacerSequenceEntry = SpacerSequence(guideEntry, metaGen)
                            data.append([spacerSequenceEntry.getGuideSequence(),
                                         spacerSequenceEntry.getOnTargetScore(),
                                         spacerSequenceEntry.getOffTargetScores(),
                                         spacerSequenceEntry.getHeuristic()])
                    if len(sys.argv) == 5:
                        saveFile = open(sys.argv[4], "w")
                        if isValidCSV(sys.argv[4]):
                            saveFile.close()
                            saveNestedListToCSV(data, sys.argv[4])
                            print(data)
                            print("Successfully saved to " + sys.argv[4])
                        elif isValidTSV(sys.argv[4]):
                            saveFile.close()
                            saveNestedListToTSV(data, sys.argv[4])
                            print(data)
                            print("Successfully saved to " + sys.argv[4])
                        else:
                            saveFile.close()
                            os.remove(sys.argv[4])
                            print(data)
                            print("[Error] Could not save data. Invalid TSV or CSV file path: " + sys.argv[4])
                    else:
                        print(data)
        else:
            print("[Error] Unknown command: " + sys.argv[1])
