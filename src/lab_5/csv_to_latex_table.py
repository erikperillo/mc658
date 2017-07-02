#!/usr/bin/env python3

import sys

def csv_to_latex_table(filepath):
    fst_line = True
    print("\\begin{table}[H]")
    print("\\centering")
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip().replace("_", "\_").split(",")
            if fst_line:
                print("\\begin{{tabular}}{{|{}|}}".format(
                    "|".join(len(line)*["c"])))
                fst_line = False
            print("\\hline")
            print(" & ".join(line) + "\\\\")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{table}")

def main():
    if len(sys.argv) < 2:
        print("usage: {} <csv_filepath>".format(sys.argv[0]))
        exit()

    csv_to_latex_table(sys.argv[1])

if __name__ == "__main__":
    main()
