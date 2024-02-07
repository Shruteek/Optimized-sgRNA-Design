import sys
import os.path
from os.path import exists
from GenomeTools import *
from MetaGenome import MetaGenome
from SpacerSequence import SpacerSequence
import time


def __fileSetup():
    """Checks if folders to be used in program exist, and if not, sets them up."""
    projectPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    outputPath = os.path.join(projectPath, "Outputs")
    dataPath = os.path.join(projectPath, "Genedata")
    if not os.path.isdir(outputPath):
        print("Project Outputs folder does not exist. Creating at " + str(outputPath))
        os.mkdir(outputPath)
    if not os.path.isdir(dataPath):
        print("Project Genedata folder does not exist. Creating at " + str(dataPath))
        os.mkdir(dataPath)


def __printHelpMessage(helpArgument):
    """Simply prints out one of the help messages explaining program details, returning nothing."""
    if helpArgument == "single-target-input" or helpArgument == "STI":
        print("Single target input format:")
        print("`python SpacerSequenceHandler.PY STI [20bp_spacer_sequence] [metagenome_file_path] [save_file_path]`\n")
        print("Single target input example:")
        print("`python SpacerSequenceHandler.PY STI TGACTGACTGACTGACTGAC metaGenome.FASTA saveFile.CSV`\n")

        print("In single target input mode, the program takes a nucleotide DNA or RNA character string representing "
              "the 20 base pair sgRNA spacer sequence of a Cas9 molecule; this string should be exactly 20 base "
              "pairs, and while it can contain 'N' nucleotides, doing so diminishes the value of the tool. The "
              "program also takes a .FASTA file path representing the metagenome environment within which the "
              "target's spacer sequence is to be analyzed, and (optionally) a .CSV or .TSV file path representing the "
              "file to which the target analysis data should be written (if the save file already exists, "
              "its contents may be overridden).\n")
        print("Each row of the target analysis data will represent an on or off-target sequence and its associated "
              "analysis data, and the columns will be as follows:")
        print("[Spacer (20 bp); "
              "\nSequence_Type (on-target or off-target); "
              "\nLocation (within the reference metagenome); "
              "\nSequence (20 bp to match the spacer); "
              "\nSurrounding_Sequence (35 bp, with the target + PAM + 6 bp before & after); "
              "\nMismatches (# of mismatches vs. spacer); "
              "\nOn_Target_Score (scaled 0 to 100); "
              "\nOff_Target_Score (scaled 0 to 100)]")
    elif helpArgument == "multi-target-input" or helpArgument == "MTI":
        print("Multi target input format:")
        print("`python SpacerSequenceHandler.PY MTI [spacers_file_path] [metagenome_file_path] [save_file_path]`\n")
        print("Multi target input example:")
        print("`python SpacerSequenceHandler.PY STI spacers.CSV metaGenome.FASTA saveFile.CSV`\n")

        print("In multi target input mode, the program takes .CSV or .TSV file path where some value represent "
              "20 base pair sgRNA spacer sequence of a Cas9 molecule; the program sifts through the file and "
              "identifies appropriate spacer sequences, and while they can contain 'N' nucleotides, "
              "doing so diminishes the value of the tool. The program also takes a .FASTA file path representing the "
              "metagenome environment within which the spacer sequences are to be analyzed, and (optionally) a "
              ".CSV or .TSV file path representing the file to which the target analysis data should be written (if "
              "the save file already exists, its contents may be overridden).\n")
        print("Each row of the target analysis data will represent an on or off-target sequence for a particular "
              "spacer, as well as its associated analysis data, and the columns will be as follows:")
        print("[Spacer (20 bp); "
              "\nSequence_Type (on-target or off-target); "
              "\nLocation (within the reference metagenome); "
              "\nSequence (20 bp to match the spacer); "
              "\nSurrounding_Sequence (35 bp, with the target + PAM + 6 bp before & after); "
              "\nMismatches (# of mismatches vs. spacer); "
              "\nOn_Target_Score (scaled 0 to 100); "
              "\nOff_Target_Score (scaled 0 to 100)]")
    else:
        print("This program can be run with any of the below input schemes:\n\n")
        print("Single target input example:")
        print("`python SpacerSequenceHandler.PY STI TGACTGACTGACTGACTGAC metaGenome.FASTA saveFile.CSV`")
        print("The program takes in a string representing a 20 base pair RNA/DNA spacer and a .FASTA metagenome file "
              "path, then saves the spacer analysis data to the saveFile path.\n")

        print("Multi target input example:")
        print("`python SpacerSequenceHandler.PY MTI targetSequences.CSV metaGenome.FASTA saveFile.CSV`")
        print("The program takes in a .CSV or .TSV file path and a .FASTA file path, sifts through the .CSV/.TSV file "
              "to identify valid 20 base pair target sequences, and saves a list of target analysis data to the "
              "saveFile path.")

        print("\n For more detailed help on specific modes, use one of the following commands:")
        print("`python SpacerSequenceHandler.py help single-target-input`")
        print("`python SpacerSequenceHandler.py help multi-target-input`""")


def __runSTI(arguments):
    """A method that runs the program given a nucleotide RNA or DNA character string representing the target
    sequence, a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to save the
    resulting data to. """
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after output file (see below): \n" + str(arguments[5:]))
    else:
        if not isValidTargetSpacerInput(arguments[2]):
            print("[Error] Invalid target sequence given: " + arguments[2] +
                  "\n (This program expects exactly 20 RNA (or DNA) base pairs representing the spacer sequence)")
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            metaGen = MetaGenome(arguments[3])
            print("Analyzing given spacer sequence: " + arguments[2])
            spacerSequence = SpacerSequence(arguments[2], metaGen)
            print("Found " + str(len(spacerSequence.getOnTargetSequences())) + " on-targets and "
                  + str(len(spacerSequence.getOffTargetSequences())) + " off-targets.")
            data = []
            appendSpacerToData(data, spacerSequence)
            print(data)
            if len(arguments) != 5: return
            projectPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
            outputPath = os.path.join(projectPath, "Outputs")
            saveFilePath = os.path.join(outputPath, os.path.basename(arguments[4]))
            saveFile = open(saveFilePath, "a+")
            if isValidCSV(saveFilePath):
                saveFile.close()
                writeNestedListToCSVRows(data, saveFilePath)
                print("Successfully saved to " + saveFilePath)
            elif isValidTSV(saveFilePath):
                saveFile.close()
                writeNestedListToTSVRows(data, saveFilePath)
                print("Successfully saved to " + saveFilePath)
            else:
                saveFile.close()
                os.remove(saveFilePath)
                print("[Error] Could not save data. Invalid TSV or CSV file path: " + saveFilePath)


def __runMTI(arguments):
    """A method that runs the program given a .CSV or .TSV file path containing character strings representing target
        sequences in cells, a .FASTA file path representing the metagenome, and (optionally) a .CSV or .TSV file path to
        save the resulting data to."""
    if len(arguments) == 2:
        print("[Error] Missing spacer sequence(s) and metagenome file inputs.")
    elif len(arguments) == 3:
        print("[Error] Missing metagenome file input.")
    elif len(arguments) > 5:
        print("[Error] Unexpected arguments after output file (see below): \n" + str(arguments[5:]))
    else:
        if not (isValidCSV(arguments[2]) or isValidTSV(arguments[2])):
            print("[Error] Invalid CSV/TSV given: " + arguments[2])
        elif not isValidFasta(arguments[3]):
            print("[Error] Invalid FASTA file path given: " + arguments[3])
        else:
            data = []
            metaGen = MetaGenome(arguments[3])
            if isValidCSV(arguments[2]):
                read_guides = csv.reader(open(arguments[2], encoding="UTF-8"), delimiter=',')
            else:
                read_guides = csv.reader(open(arguments[2], encoding="UTF-8"), delimiter="\t")
            for row in read_guides:
                for entry in row:
                    if isValidTargetSpacerInput(entry):
                        print("Analyzing given spacer sequence: " + entry)
                        spacerSequenceEntry = SpacerSequence(entry, metaGen)
                        print("Found " + str(len(spacerSequenceEntry.getOnTargetSequences())) + " on-targets and "
                              + str(len(spacerSequenceEntry.getOffTargetSequences())) + " off-targets.")
                        appendSpacerToData(data, spacerSequenceEntry)
            print(data)
            if len(arguments) != 5: return
            saveFilePath = os.path.join(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                                     "Outputs"), os.path.basename(arguments[4]))
            saveFile = open(saveFilePath, "a+")
            saveFile.close()
            if isValidCSV(saveFilePath) or isValidTSV(saveFilePath):
                writeNestedListToRows(data, saveFilePath)
                print("Successfully saved to " + saveFilePath)
            else:
                os.remove(saveFilePath)
                print("[Error] Could not save data. Invalid TSV or CSV file path: " + saveFilePath)


if __name__ == '__main__':
    """Handles SeedSequence.py, accepting inputs and outputs to the class and returning the result of on-target and 
    off-target calculations and heuristics."""
    __fileSetup()
    if len(sys.argv) == 1:
        __printHelpMessage("")
    elif len(sys.argv) >= 2:
        startTime = time.time()
        if sys.argv[1] == "help":
            if len(sys.argv) == 2:
                __printHelpMessage("")
            elif len(sys.argv) == 3:
                __printHelpMessage(sys.argv[2])
            elif len(sys.argv) > 3:
                print("[Error] Unexpected arguments after help argument: " + str(sys.argv[3:]))
        elif sys.argv[1] == "STI":
            __runSTI(sys.argv)
        elif sys.argv[1] == "MTI":
            __runMTI(sys.argv)
        else:
            print("[Error] Unknown command: " + sys.argv[1])
        print("Program runtime (seconds): " + str(time.time() - startTime))

