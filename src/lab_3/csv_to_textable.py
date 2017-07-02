#!/usr/bin/env python3

import sys

def main():
    if len(sys.argv) < 2:
        return

    with open(sys.argv[1]) as f:
        for l in f:
            l = l.strip()
            l = l.replace(",", " & ")
            l = l.replace("_", "\_")
            print(l, "\\\\", sep="")
            print("\hline")

if __name__ == "__main__":
    main()
