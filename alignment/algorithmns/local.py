import sys
import pandas as pd
import numpy as np
from utilities import DIAGONAL, UP, LEFT, generate_scoring_matrix


def initialize_alignment_table(file):
    """Initializes the alignment table and pointer table for sequence alignment.
    Args:
        file (str): The path to a file containing two sequences to be aligned, each on a separate line.
    Returns:
        tuple: A tuple containing the initialized alignment table (DataFrame), pointer table (DataFrame),
               the first sequence (str), and the second sequence (str).  The alignment table and pointer
               table are Pandas DataFrames. The alignment table is initialized with empty values,
               and the first row and column are filled with zeros. The pointer table is initialized
               with zeros.
    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file does not contain exactly two lines.
    """
    try:
        with open(file, "r") as text:
            text = text.readlines()
    except FileNotFoundError:
        print(f"File {file} not found.")
        sys.exit(1)

    try:
        if len(text) != 2:
            raise ValueError("File must contain two lines.")
    except ValueError as e:
        print(e)
        sys.exit(1)

    first_sequence = text[0].strip()
    second_sequence = text[1].strip()
    first_line = ["-"] + list(first_sequence)
    second_line = ["-"] + list(second_sequence)

    alignment_table = pd.DataFrame(np.empty(
        (len(second_line), len(first_line))))
    alignment_table.columns = first_line
    alignment_table.index = second_line

    first_row_values = np.zeros((len(first_sequence)), dtype=int)
    second_row_values = np.zeros((len(second_sequence)), dtype=int)
    alignment_table.iloc[0, 1:] = first_row_values
    alignment_table.iloc[1:, 0] = second_row_values

    pointer_table = pd.DataFrame(0,
                                 index=alignment_table.index,
                                 columns=alignment_table.columns)
    pointer_table.iloc[0, 1:] = 0
    pointer_table.iloc[1:, 0] = 0

    return alignment_table, pointer_table, first_sequence, second_sequence


def smith_waterman(alignment_file, scoring_file, gap_penalty):
    """
    Performs local sequence alignment using the Smith-Waterman algorithm.
    Args:
        alignment_file (str): Path to the file containing the sequences to align.
        scoring_file (str): Path to the file containing the scoring matrix for matches and mismatches.
        gap_penalty (int): The penalty for introducing a gap in the alignment.
    Returns:
        tuple: A tuple containing the maximum alignment score, the aligned sequence 1, and the aligned sequence 2.
               (max_score, aligned_seq1, aligned_seq2)
    """

    scoring_table = generate_scoring_matrix(scoring_file)
    alignment_table, pointer_table, first_sequence, second_sequence = (
        initialize_alignment_table(alignment_file))
    for index in range(1, len(alignment_table.index)):
        for column in range(1, len(alignment_table.columns)):
            match = int(scoring_table[alignment_table.index[index]][
                alignment_table.columns[column]])
            diagonal = int(alignment_table.iloc[index - 1, column - 1]) + match
            up = int(alignment_table.iloc[index - 1, column]) + gap_penalty
            left = int(alignment_table.iloc[index, column - 1]) + gap_penalty
            max_score = max(diagonal, up, left, 0)
            alignment_table.iloc[index, column] = max_score

            direction = 0
            if max_score == 0:
                direction = 0
            else:
                if diagonal == max_score:
                    direction |= DIAGONAL
                if up == max_score:
                    direction |= UP
                if left == max_score:
                    direction |= LEFT

            pointer_table.iloc[index, column] = direction
    max_score, rows, columns = find_max_score(alignment_table)
    aligned_seq1, aligned_seq2 = backtracking(first_sequence, second_sequence,
                                              alignment_table, pointer_table)

    return max_score, aligned_seq1, aligned_seq2


def backtracking(seq1, seq2, alignment_table, pointer_table):
    """
    Perform backtracking to reconstruct the optimal local alignment from the alignment and pointer tables.
    The function starts from the cell with the highest score in the alignment table and traces back
    using the pointer table to reconstruct the aligned sequences. It stops when a cell with a pointer
        value of 0 is reached, indicating the start of the alignment.
        Args:
                seq1 (str): The first sequence.
                seq2 (str): The second sequence.
                alignment_table (pandas.DataFrame): The table containing alignment scores.
                pointer_table (pandas.DataFrame): The table containing pointers indicating the direction of movement.
        Returns:
                tuple: A tuple containing the aligned sequences (aligned_seq1, aligned_seq2).
                           The aligned sequences are built by tracing back from the cell with the highest score
                           to the beginning of the alignment based on the pointer values.
        Raises:
        ValueError: If an invalid pointer value is encountered during backtracking.
    """

    # Start at cell with highet value --> need to get column and row values
    max_score, rows, columns = find_max_score(alignment_table)

    aligned_seq1 = ""
    aligned_seq2 = ""

    current_pointer = pointer_table.iloc[rows, columns]

    while current_pointer != 0:
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
        current_pointer = pointer_table.iloc[rows, columns]

    return aligned_seq1, aligned_seq2
    # return aligned_score, aligned_seq1, aligned_seq2


def find_max_score(alignment_table):
    """Finds the maximum score in the alignment table.

    Args:
            alignment_table (pandas.DataFrame): The alignment table.

    Returns:
            tuple: A tuple containing the maximum score and its position (row, column).
    """
    max_score = alignment_table.max().max()
    row, column = np.unravel_index(alignment_table.values.argmax(),
                                   alignment_table.shape)

    return max_score, row, column


if __name__ == "__main__":

    # Example usage
    alignment_file = "alignment/test/sequences/03.txt"
    scoring_file = "alignment/test/matrices/standard.m"
    gap_penalty = -1
    print(smith_waterman(alignment_file, scoring_file, gap_penalty))
