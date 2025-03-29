import sys
import pandas as pd
import numpy as np


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
        first_line = list(text[0].strip())
        first_line.insert(0, "-")
        second_line = list(text[1].strip())
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

    return alignment_table


def needleman_wunsch(alignment_file, scoring_file, gap_penalty):
    scoring_table = build_scoring_matrix(scoring_file)
    alignment_table = initialize_alignment_table(alignment_file, gap_penalty)

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
            alignment_table.iloc[index, column] = max(diagonal, up, left)
    print(alignment_table)


needleman_wunsch("01.txt", "standard.m", -1)



