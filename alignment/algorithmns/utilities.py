import sys 
import numpy as np
import pandas as pd


# Constants for pointer directions
DIAGONAL = 1  # 001
UP = 2        # 010
LEFT = 4      # 100


def generate_scoring_matrix(file):
    """Builds a scoring matrix from a file.

    The file should be a text file with the scoring matrix in the following format:
    The first line should contain the nucleotides, separated by spaces.
    The following lines should contain the scores for each nucleotide, with the first element being the nucleotide and the following elements being the scores for each nucleotide in the order specified in the first line.
    
    Args:
		file (str): The path to the file containing the scoring matrix.
	Raises:
		FileNotFoundError: If the file does not exist.
		ValueError: If the file is not in the correct format.
	Returns:
		dict: A dictionary representing the scoring matrix. The keys are the nucleotides, and the values are dictionaries containing the scores for each nucleotide.
    """
    scoring_dict = {}
    try:
        with open(file, "r") as text:
            text = text.readlines()
    except FileNotFoundError:
        print(f"File {file} not found.")
        sys.exit(1)
    
    # TODO: Add error handling for incorrect file format
		# Check if the file has at least two lines
		# Check matrix dimensions - number of rows and columns
		# Check that first line has correct header placeholder as first element 
		# Check each row element matches number of columns 
    
    nucleotides = text[0].strip().split(" ")[1:]

    for lines in text[1:]:
        lines = lines.strip().split()
        scoring_dict[lines[0]] = {
            nucleotides[index - 1]: lines[index]
            for index in range(1,
                               len(nucleotides) + 1)
        }

    return scoring_dict
