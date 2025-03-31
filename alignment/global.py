import sys
import pandas as pd
import numpy as np
from alignment_constants import DIAGONAL, UP, LEFT


def build_scoring_matrix(file):
    scoring_dict = {}
    with open(file, "r") as text:
        text = text.readlines()
        nucleotides = text[0].strip().split(" ")[1:]
    for lines in text[1:]:
        lines = lines.strip().split()
        scoring_dict[lines[0]] = {
            nucleotides[index - 1]: lines[index]
            for index in range(1, len(nucleotides) + 1)
        }
    return scoring_dict


def initialize_alignment_table(file, gap_penalty):
    with open(file, "r") as text:
        text = text.readlines()
        first_sequence = text[0].strip()
        second_sequence = text[1].strip()
        first_line = list(first_sequence)
        first_line.insert(0, "-")
        second_line = list(second_sequence)
        second_line.insert(0, "-")
    alignment_table = pd.DataFrame(np.empty((len(second_line), len(first_line))))
    alignment_table.columns = first_line
    alignment_table.index = second_line
    first_row_values = np.arange(
        gap_penalty, (len(first_line)) * gap_penalty, gap_penalty
    )
    second_row_values = np.arange(
        gap_penalty, (len(second_line)) * gap_penalty, gap_penalty
    )
    alignment_table.iloc[0, 1:] = first_row_values
    alignment_table.iloc[1:, 0] = second_row_values

    pointer_table = pd.DataFrame(
        0, index=alignment_table.index, columns=alignment_table.columns
    )
    pointer_table.iloc[0, 1:] = LEFT
    pointer_table.iloc[1:, 0] = UP

    return alignment_table, pointer_table, first_sequence, second_sequence


def needleman_wunsch(alignment_file, scoring_file, gap_penalty):
    scoring_table = build_scoring_matrix(scoring_file)
    alignment_table, pointer_table, first_sequence, second_sequence = (
        initialize_alignment_table(alignment_file, gap_penalty)
    )

    for index in range(1, len(alignment_table.index)):
        for column in range(1, len(alignment_table.columns)):
            match = int(
                scoring_table[alignment_table.index[index]][
                    alignment_table.columns[column]
                ]
            )
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
    aligned_seq1, aligned_seq2 = backtracking(
        first_sequence, second_sequence, pointer_table
    )
    aligned_score = alignment_table.iloc[-1, -1]
    return aligned_score, aligned_seq1, aligned_seq2
   


def backtracking(seq1, seq2, pointer_table):
    # Start at bottom-right cell
    columns = len(pointer_table.columns) - 1 
    rows = len(pointer_table.index) - 1 

    aligned_seq1 = ""
    aligned_seq2 = ""

    while rows > 0 or columns > 0:
        current_pointer = pointer_table.iloc[rows,columns]

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
