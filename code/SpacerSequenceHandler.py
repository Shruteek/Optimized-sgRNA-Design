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


def __queryForGuideSequence():
    """Repeatedly queries the user for a string guide sequence, validates each sequence input, takes the first valid
    sequence, and returns it."""
    guideSequence = ""
    timesQueried = 0
    while not (
            isValidRNA(guideSequence) and (len(guideSequence) == 20 or len(guideSequence) == 35)) and timesQueried < 5:
        if not timesQueried == 0:
            print("Invalid guide sequence or spacer sequence entered. Please try again.")
        guideSequence = input("Enter a 20 bp spacer sequence: ")
        timesQueried = timesQueried + 1
        if timesQueried == 5 and \
                not (isValidRNA(guideSequence) and (len(guideSequence) == 20 or len(guideSequence) == 35)):
            print("Queried too many times. Exiting...")
    if not (isValidRNA(guideSequence) and (len(guideSequence) == 20 or len(guideSequence) == 35)):
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
            print("Invalid save file path entered. Please try again.")
        filePathInput = input("Enter the desired save file path (e.g. saveFile.CSV or /c/Users/saveFile.TSV), "
                              "including the type extension: ")
        if len(filePathInput) > 4 and (str.lower(filePathInput[-4:]) == ".tsv" or
                                       str.lower(filePathInput[-4:]) == ".csv"):
            validInput = True
        timesQueried = timesQueried + 1
        if timesQueried == 5 and not validInput:
            print("Queried too many times. Exiting...")
    if not validInput:
        return ""
    else:
        print(filePathInput)
        return filePathInput


def __printHelpMessage(helpArgument):
    """Simply prints out one of the help messages explaining program details, returning nothing."""
    if helpArgument == "no-input":
        print("""No input (e.g. `python SpacerSequenceHandler.py`):\nIn no input mode, the program can take a 20 
        nucleotide character string representing the 20 base pair sgRNA spacer sequence and a metagenome .FASTA file 
        and return the target analysis list of the sgRNA. It can also take a delimiter-type file (e.g. CSV or TSV) of 
        sgRNA spacers and a metagenome .FASTA file and return a list of target analysis lists, one for each valid 
        sgRNA. These inputs are repeatedly queried for command line text input until 5 invalid inputs are given for a 
        prompt, or until all valid inputs have been given.\nNOTICE: The no-input mode is deprecated and its entire 
        functionality can be executed through single command line commands.""")
    elif helpArgument == "single-spacer-input" or helpArgument == "SSI":
        print("""Single spacer input (e.g `python SpacerSequenceHandler.PY SSI ACUGACUGACUGACUGACUG metaGenome.FASTA 
        saveFile.CSV`):\nIn single spacer input mode, the program can take a nucleotide RNA or DNA character string 
        representing the 20 base pair sgRNA spacer sequence; this string can be 20 base pairs (just the spacer), 
        23 base pairs (the spacer + PAM complement), or 35 base pairs (the guide sequence containing the spacer), 
        and in RNA or DNA format (though currently, no generic 'N' functionality). The program also takes a .FASTA 
        file path representing the metaGenome with which the spacer sequence is to be analyzed, and (optionally) a 
        .CSV or .TSV file path representing the file to which the target analysis data should be written (**if the 
        save file already exists, its contents will not be appended to**). The target analysis data will take the 
        form:\n[ Spacer (20 bp RNA), Target Guide Sequence (35 bp DNA), Off-Target Sequences (35 bp DNA), On-Target Score, 
        Off-Target Scores, Heuristic]""")
    elif helpArgument == "single-target-input" or helpArgument == "STI":
        print("""Single target input (e.g `python SpacerSequenceHandler.PY STI TGACTGACTGACTGACTGAC metaGenome.FASTA 
        saveFile.CSV`):\nIn single target input mode, the program can take a nucleotide DNA character string 
        representing the 20 base pair sgRNA spacer sequence; this string can be 20 base pairs (just the target), 
        23 base pairs (the target + PAM), or 35 base pairs (the guide target sequence containing the target), 
        only in DNA format (though currently, no generic 'N' functionality). The program also takes a .FASTA file 
        path representing the metaGenome with which the target's spacer sequence is to be analyzed, and (optionally) 
        a .CSV or .TSV file path representing the file to which the target analysis data should be written (**if the 
        save file already exists, its contents will not be appended to**). The target analysis data will take the 
        form:\n[ Spacer (20 bp RNA), Target Guide Sequence (35 bp DNA), Off-Target Sequences (35 bp DNA), On-Target Score, 
        Off-Target Scores, Heuristic]""")
    elif helpArgument == "multi-spacer-input" or helpArgument == "MSI":
        print("""Multi spacer input (e.g `python SpacerSequenceHandler.PY MSI spacerSequences.CSV metaGenome.FASTA 
        saveFile.CSV`):\nIn multi spacer input mode, the program can take a .TSV or .CSV file path representing the 
        file from which the sgRNA string sequences will be read; these strings can be 20 base pairs (just the 
        spacers), 23 base pairs (the spacers + NGG PAM complements), or 35 base pairs (the guide sequences containing 
        the spacers), and in RNA or DNA format (though currently, no generic 'N' functionality). The program also 
        takes in a .FASTA file path representing the metaGenome with which the spacer sequences are to be analyzed, 
        and (optionally) a .CSV or .TSV file path representing the file to which the target analysis data should be 
        written (**if the save file already exists, its contents will not be appended to**). Each target analysis row 
        will take the form:\n[Spacer (20 bp RNA), Target Guide Sequence (35 bp DNA), Off-Target Sequences (35 bp DNA), 
        On-Target Score, Off-Target Scores, Heuristic]""")
    elif helpArgument == "multi-target-input" or helpArgument == "MTI":
        print("""Multi target input (e.g `python SpacerSequenceHandler.PY MTI targetSequences.CSV metaGenome.FASTA 
        saveFile.CSV`):\nIn multi target input mode, the program can take a .TSV or .CSV file path representing the 
        file from which the DNA string target sequences will be read; these strings can be 20 base pairs (just the 
        targets), 23 base pairs (the targets + NGG PAMs), or 35 base pairs (the target guide sequences containing the 
        target), and in RNA or DNA format (though currently, no generic 'N' functionality). The program also takes in 
        a .FASTA file path representing the metaGenome with which the target sequences are to be found and analyzed, 
        and (optionally) a .CSV or .TSV file path representing the file to which the target analysis data should be 
        written (**if the save file already exists, its contents will not be appended to**). Each target analysis row 
        will take the form:\n[Spacer (20 bp RNA), Target Guide Sequence (35 bp DNA), Off-Target Sequences (35 bp DNA), 
        On-Target Score, Off-Target Scores, Heuristic]""")
    else:
        if not helpArgument == "":
            print("[Error] Invalid help argument Entered. See generic help message below.")
        print("""This program can be run with any of the below input schemes:\n\nNo input (e.g. `python 
        SpacerSequenceHandler.PY`): The program will repeatedly prompt for inputs until it receives all necessary 
        information or too many failed inputs.\n\nSingle spacer input (e.g `python SpacerSequenceHandler.PY SSI 
        ACUGACUGACUGACUGACUG metaGenome.FASTA saveFile.CSV`): The program takes in a string representing an sgRNA 
        spacer and a .FASTA file path and saves the target analysis data to the saveFile path.\n\nSingle target input 
        (e.g `python SpacerSequenceHandler.PY STI TGACTGACTGACTGACTGAC metaGenome.FASTA saveFile.CSV`): The program 
        takes in a string representing a DNA target and a .FASTA file path and saves the target analysis data to the 
        saveFile path.\n\nMulti spacer input (e.g `python SpacerSequenceHandler.PY MSI spacerSequences.CSV 
        metaGenome.FASTA saveFile.CSV`): The program takes in a .CSV or .TSV file path containing spacer sequences 
        and a .FASTA file path and saves a list of target analysis data to the saveFile path.\n\nMulti target input ( 
        e.g `python SpacerSequenceHandler.PY MTI targetSequences.CSV metaGenome.FASTA saveFile.CSV`):The program 
        takes in a .CSV or .TSV file path containing target sequences and a .FASTA file path and saves a list of 
        target analysis data to the saveFile path.\n\nHelp (e.g. `python SpacerSequenceHandler.py help`): Print this 
        message. For more detailed help on specific modes, use one of the following commands:\n`python 
        SpacerSequenceHandler.py help no-input`\n`python SpacerSequenceHandler.py help single-spacer-input`\n`python 
        SpacerSequenceHandler.py help single-target-input`\n`python SpacerSequenceHandler.py help 
        multi-spacer-input`\n`python SpacerSequenceHandler.py help multi-target-input`""")


def __runNoInput():
    """Runs the program, from deciding the program mode to returning an output."""
    print("You have entered no input mode. This program can be run with other input")
    inputMode = __queryForProgramType()
    if inputMode == 1 or inputMode == 2:
        metaGenPath = __queryForMetagenomeFile()
        metaGenome = MetaGenome(metaGenPath, metaGenPath)
        headerRow = ["Spacer (20 bp RNA)", "Target Guide Sequence (35 bp DNA)", "Guide Sequence Heuristic",
                     "On-Target Score", "Off-Targets (35 bp DNA)", "Off-Target Scores", "Off-Target Count"]
        data = []
        if metaGenome.size() == 0:
            print("Empty or invalid metagenome provided. Exiting...")
            return
        elif inputMode == 1:
            guideSequence = __queryForGuideSequence()
            if not (len(guideSequence) == 20 or isValidGuideSequence(guideSequence)):
                print("No valid guide sequences (or contained spacer sequences) given. Exiting...")
                return
            spacerSequence = SpacerSequence(guideSequence, metaGenome)
            for guideSequenceIndex in range(len(spacerSequence.getGuideSequences())):
                data.append([spacerSequence.getSpacerSequence(),
                             complementaryDNA(spacerSequence.getGuideSequences()[guideSequenceIndex]),
                             spacerSequence.getOffTargetSequences(),
                             spacerSequence.getOnTargetScores()[guideSequenceIndex],
                             spacerSequence.getOffTargetScores()[guideSequenceIndex],
                             spacerSequence.getHeuristics()[guideSequenceIndex]])

        elif inputMode == 2:
            delimiter = __queryForDelimiterType()
            delimiterFile = __queryForDelimitedGuideFile()
            if len(delimiterFile) == 0:
                print("No valid delimiter files given. Exiting...")
                return
            guidesTable = csv.reader(open(delimiterFile, newline=''), delimiter=delimiter)
            for guideRow in guidesTable:
                for guideColumnEntry in guideRow:
                    spacerSequenceEntry = SpacerSequence(guideColumnEntry, metaGenome)
                    for guideSequenceIndex in range(len(spacerSequenceEntry.getGuideSequences())):
                        data.append([spacerSequenceEntry.getSpacerSequence(),
                                     complementaryDNA(spacerSequenceEntry.getGuideSequences()[guideSequenceIndex]),
                                     spacerSequenceEntry.getOffTargetSequences(),
                                     spacerSequenceEntry.getOnTargetScores()[guideSequenceIndex],
                                     spacerSequenceEntry.getOffTargetScores()[guideSequenceIndex],
                                     spacerSequenceEntry.getHeuristics()[guideSequenceIndex]])

        doSaveToFile = __queryForSave()
        if doSaveToFile:
            saveFilePath = __queryForSaveFile()
            if not len(saveFilePath) == 0:
                saveFile = open(saveFilePath, "a+")
                saveFile.close()
                if isValidCSV(saveFilePath):
                    writeListToCSVRow(headerRow, saveFilePath)
                    writeNestedListToCSVRows(data, saveFilePath)
                    print("Program data saved to " + saveFilePath)
                elif isValidTSV(saveFilePath):
                    writeListToTSVRow(headerRow, saveFilePath)
                    writeNestedListToTSVRows(data, saveFilePath)
                    print("Program data saved to " + saveFilePath)
    elif inputMode == 3:
        print("Exiting program...")
        return


def __runSSI(arguments):
    """A method that runs the program given a nucleotide RNA or DNA character string representing the spacer sequence,
    a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to save the resulting data
    to."""
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after metagenome file input: " + str(arguments[5:]))
    else:
        if not isValidSpacerInput(arguments[2]):
            print("[Error] Invalid sgRNA sequence given: " + arguments[2] +
                  " (expected 20, 23, or 35 DNA/RNA base pairs)")
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            headerRow = ["Spacer (20 bp RNA)", "Target Guide Sequence (35 bp DNA)", "Guide Sequence Heuristic",
                         "On-Target Score", "Off-Targets (35 bp DNA)", "Off-Target Scores", "Off-Target Count"]
            data = []
            metaGen = MetaGenome(arguments[3])
            print("Analyzing given spacer sequence: " + complementaryRNA(complementaryRNA(arguments[2])))
            spacerSequence = SpacerSequence(complementaryRNA(complementaryRNA(arguments[2])), metaGen)
            appendSpacerToData(data, spacerSequence)
            print(data)
            if len(arguments) == 5:
                saveFile = open(arguments[4], "a+")
                if isValidCSV(arguments[4]):
                    saveFile.close()
                    writeListToCSVRow(headerRow, arguments[4])
                    writeNestedListToCSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                elif isValidTSV(arguments[4]):
                    saveFile.close()
                    writeListToTSVRow(headerRow, arguments[4])
                    writeNestedListToTSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                else:
                    saveFile.close()
                    os.remove(arguments[4])
                    print("[Error] Could not save data. Invalid TSV or CSV file path: " + arguments[4])


def __runSTI(arguments):
    """A method that runs the program given a nucleotide RNA or DNA character string representing the target
    sequence, a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to save the
    resulting data to. """
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after metagenome file input: " + str(arguments[5:]))
    else:
        if not isValidTargetInput(arguments[2]):
            print("[Error] Invalid target sequence given: " + arguments[2] +
                  " (expected 20, 23, or 35 DNA base pairs)")
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            metaGen = MetaGenome(arguments[3])
            print("Analyzing given spacer sequence: " + complementaryRNA(arguments[2]))
            spacerSequence = SpacerSequence(complementaryRNA(arguments[2]), metaGen)
            headerRow = ["Spacer (20 bp RNA)", "Target Guide Sequence (35 bp DNA)", "Guide Sequence Heuristic",
                         "On-Target Score", "Off-Targets (35 bp DNA)", "Off-Target Scores", "Off-Target Count"]
            data = []
            appendSpacerToData(data, spacerSequence)
            print(data)
            if len(arguments) == 5:
                saveFile = open(arguments[4], "a+")
                if isValidCSV(arguments[4]):
                    saveFile.close()
                    writeListToCSVRow(headerRow, arguments[4])
                    writeNestedListToCSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                elif isValidTSV(arguments[4]):
                    saveFile.close()
                    writeListToTSVRow(headerRow, arguments[4])
                    writeNestedListToTSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                else:
                    saveFile.close()
                    os.remove(arguments[4])
                    print("[Error] Could not save data. Invalid TSV or CSV file path: " + arguments[4])


def __runMSI(arguments):
    """A method that runs the program given a .CSV or .TSV file path containing character strings representing spacer
    sequences in cells, a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to
    save the resulting data to."""
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after metagenome file input: " + str(arguments[5:]))
    else:
        if not (isValidCSV(arguments[2]) or isValidTSV(arguments[2])):
            print("[Error] Invalid CSV/TSV given: " + arguments[2])
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            data = []
            headerRow = ["Spacer (20 bp RNA)", "Target Guide Sequence (35 bp DNA)", "Guide Sequence Heuristic",
                         "On-Target Score", "Off-Targets (35 bp DNA)", "Off-Target Scores", "Off-Target Count"]
            metaGen = MetaGenome(arguments[3])
            if isValidCSV(arguments[2]):
                read_guides = csv.reader(open(arguments[2], newline=''), delimiter=',')
            else:
                read_guides = csv.reader(open(arguments[2], newline=''), delimiter="\t")
            for row in read_guides:
                for entry in row:
                    if isValidSpacerInput(entry):
                        print("Analyzing given spacer sequence: " + complementaryRNA(complementaryRNA(entry)))
                        spacerSequenceEntry = SpacerSequence(complementaryRNA(complementaryRNA(entry)), metaGen)
                        appendSpacerToData(data, spacerSequenceEntry)
            print(data)
            if len(arguments) == 5:
                saveFile = open(arguments[4], "a+")
                if isValidCSV(arguments[4]):
                    saveFile.close()
                    writeListToCSVRow(headerRow, arguments[4])
                    writeNestedListToCSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                elif isValidTSV(arguments[4]):
                    saveFile.close()
                    writeListToTSVRow(headerRow, arguments[4])
                    writeNestedListToTSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                else:
                    saveFile.close()
                    os.remove(arguments[4])
                    print("[Error] Could not save data. Invalid TSV or CSV file path: " + arguments[4])


def __runMTI(arguments):
    """A method that runs the program given a .CSV or .TSV file path containing character strings representing target
        sequences in cells, a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to
        save the resulting data to."""
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after metagenome file input: " + str(arguments[5:]))
    else:
        if not (isValidCSV(arguments[2]) or isValidTSV(arguments[2])):
            print("[Error] Invalid CSV/TSV given: " + arguments[2])
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            data = []
            headerRow = ["Spacer (20 bp RNA)", "Target Guide Sequence (35 bp DNA)", "Guide Sequence Heuristic",
                         "On-Target Score", "Off-Targets (35 bp DNA)", "Off-Target Scores", "Off-Target Count"]
            metaGen = MetaGenome(arguments[3])
            if isValidCSV(arguments[2]):
                read_guides = csv.reader(open(arguments[2], newline=''), delimiter=',')
            else:
                read_guides = csv.reader(open(arguments[2], newline=''), delimiter="\t")
            for row in read_guides:
                for entry in row:
                    if isValidTargetInput(entry):
                        print("Analyzing given spacer sequence: " + complementaryRNA(entry))
                        spacerSequenceEntry = SpacerSequence(complementaryRNA(entry), metaGen)
                        appendSpacerToData(data, spacerSequenceEntry)
            print(data)
            if len(arguments) == 5:
                saveFile = open(arguments[4], "a+")
                if isValidCSV(arguments[4]):
                    saveFile.close()
                    writeListToCSVRow(headerRow, arguments[4])
                    writeNestedListToCSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                elif isValidTSV(arguments[4]):
                    saveFile.close()
                    writeListToTSVRow(headerRow, arguments[4])
                    writeNestedListToTSVRows(data, arguments[4])
                    print("Successfully saved to " + arguments[4])
                else:
                    saveFile.close()
                    os.remove(arguments[4])
                    print("[Error] Could not save data. Invalid TSV or CSV file path: " + arguments[4])


if __name__ == '__main__':
    """Handles SeedSequence.py, accepting inputs and outputs to the class and returning the result of on-target and 
    off-target calculations and heuristics."""
    if len(sys.argv) == 1:
        __runNoInput()
    elif len(sys.argv) >= 2:
        if sys.argv[1] == "help":
            if len(sys.argv) == 2:
                __printHelpMessage("")
            elif len(sys.argv) == 3:
                __printHelpMessage(sys.argv[2])
            elif len(sys.argv) > 3:
                print("[Error] Unexpected arguments after help argument: " + str(sys.argv[3:]))
        elif sys.argv[1] == "SSI":
            __runSSI(sys.argv)
        elif sys.argv[1] == "STI":
            __runSTI(sys.argv)
        elif sys.argv[1] == "MSI":
            __runMSI(sys.argv)
        elif sys.argv[1] == "MTI":
            __runMTI(sys.argv)
        else:
            print("[Error] Unknown command: " + sys.argv[1])
