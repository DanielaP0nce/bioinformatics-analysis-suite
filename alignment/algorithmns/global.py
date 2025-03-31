import sys
import pandas as pd
import numpy as np
from alignment.algorithmns.alignment_utilities import DIAGONAL, UP, LEFT


def build_scoring_matrix(file):
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
    nucleotides = text[0].strip().split(" ")[1:]

    for lines in text[1:]:
        lines = lines.strip().split()
        scoring_dict[lines[0]] = {
            nucleotides[index - 1]: lines[index]
            for index in range(1,
                               len(nucleotides) + 1)
        }

    return scoring_dict


def initialize_alignment_table(file, gap_penalty):
    """Initializes the alignment table and pointer table for the Needleman-Wunsch algorithm.
        Args:
            file (str): The path to the file containing the two sequences to be aligned. The file should contain two lines, each representing a sequence.
            gap_penalty (int): The penalty for introducing a gap in the alignment.
        Returns:
            tuple: A tuple containing the initialized alignment table (pandas.DataFrame), the pointer table (pandas.DataFrame), the first sequence (str), and the second sequence (str).  The alignment table is initialized with gap penalties along the first row and first column. The pointer table is initialized with 'LEFT' for the first row (excluding the first element) and 'UP' for the first column (excluding the first element).  All other elements are initialized to 0.
        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file does not contain two lines.
    """
    try:
        with open(file, "r") as text:
            text = text.readlines()
    except FileNotFoundError:
        print(f"File {file} not found.")
        sys.exit(1)

    first_sequence = text[0].strip()
    second_sequence = text[1].strip()

    first_line = list(first_sequence)
    first_line.insert(0, "-")
    second_line = list(second_sequence)
    second_line.insert(0, "-")

    alignment_table = pd.DataFrame(np.empty(
        (len(second_line), len(first_line))))
    alignment_table.columns = first_line
    alignment_table.index = second_line

    first_row_values = np.arange(gap_penalty, (len(first_line)) * gap_penalty,
                                 gap_penalty)
    second_row_values = np.arange(gap_penalty, (len(second_line)) * gap_penalty,
                                  gap_penalty)
    alignment_table.iloc[0, 1:] = first_row_values
    alignment_table.iloc[1:, 0] = second_row_values

    pointer_table = pd.DataFrame(0,
                                 index=alignment_table.index,
                                 columns=alignment_table.columns)
    pointer_table.iloc[0, 1:] = LEFT
    pointer_table.iloc[1:, 0] = UP

    return alignment_table, pointer_table, first_sequence, second_sequence


def needleman_wunsch(alignment_file, scoring_file, gap_penalty):
    """
     Performs global alignment of two sequences using the Needleman-Wunsch algorithm.
     Args:
      alignment_file (str): Path to the file containing the sequences to be aligned.
      scoring_file (str): Path to the file containing the scoring matrix for matches and mismatches.
      gap_penalty (int): The penalty for introducing a gap in the alignment.
     Returns:
      tuple: A tuple containing the alignment score, the aligned first sequence, and the aligned second sequence.
          The alignment score is an integer representing the optimal score found by the algorithm.
          The aligned sequences are strings representing the sequences after introducing gaps for optimal alignment.
     Raises:
      FileNotFoundError: If the specified files do not exist."""
    # Check if the files exist

    scoring_table = build_scoring_matrix(scoring_file)
    alignment_table, pointer_table, first_sequence, second_sequence = (
        initialize_alignment_table(alignment_file, gap_penalty))

    for index in range(1, len(alignment_table.index)):
        for column in range(1, len(alignment_table.columns)):
            match = int(scoring_table[alignment_table.index[index]][
                alignment_table.columns[column]])
            diagonal = int(alignment_table.iloc[index - 1, column - 1]) + match
            up = int(alignment_table.iloc[index - 1, column]) + gap_penalty
            left = int(alignment_table.iloc[index, column - 1]) + gap_penalty
            max_score = max(diagonal, up, left)
            alignment_table.iloc[index, column] = max_score

            direction = 0
            if diagonal == max_score:
                direction |= DIAGONAL
            if up == max_score:
                direction |= UP
            if left == max_score:
                direction |= LEFT

            pointer_table.iloc[index, column] = direction

    aligned_seq1, aligned_seq2 = backtracking(first_sequence, second_sequence,
                                              pointer_table)
    aligned_score = alignment_table.iloc[-1, -1]

    return aligned_score, aligned_seq1, aligned_seq2


def backtracking(seq1, seq2, pointer_table):
    """Performs backtracking on a pointer table to reconstruct the aligned sequences.
        The function starts from the bottom-right cell of the pointer table and
        traces back the optimal alignment path based on the pointers stored in each cell.
        It reconstructs the aligned sequences by following the path and inserting gaps ('-')
        where necessary.
		Args:
			seq1 (str): The first sequence.
			seq2 (str): The second sequence.
			pointer_table (pandas.DataFrame): A DataFrame representing the pointer table,
				where each cell contains a bitmask indicating the direction(s) from which
				the optimal alignment score was derived (DIAGONAL, UP, LEFT).
		Returns:
			tuple: A tuple containing two strings, representing the aligned versions
				of seq1 and seq2, respectively.
		Raises:
            ValueError: If an invalid pointer value is encountered during backtracking
                (i.e., a cell contains a pointer value that does not correspond to
                a valid direction).
        """

    # Start at bottom-right cell
    columns = len(pointer_table.columns) - 1
    rows = len(pointer_table.index) - 1

    aligned_seq1 = ""
    aligned_seq2 = ""

    while rows > 0 or columns > 0:
        current_pointer = pointer_table.iloc[rows, columns]

        if current_pointer & DIAGONAL and rows > 0 and columns > 0:
            # Diagonal move (match/mismatch)
            aligned_seq1 = seq1[columns - 1] + aligned_seq1
            aligned_seq2 = seq2[rows - 1] + aligned_seq2
            rows -= 1
            columns -= 1
        elif current_pointer & UP and rows > 0:
            # Up move (gap in seq1)
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[rows - 1] + aligned_seq2
            rows -= 1
        elif current_pointer & LEFT and columns > 0:
            # Up move (gap in seq2)
            aligned_seq1 = seq1[columns - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            columns -= 1
        else:
            raise ValueError("Invalid pointer value")

    return aligned_seq1, aligned_seq2


needleman_wunsch("01.txt", "standard.m", -1)
