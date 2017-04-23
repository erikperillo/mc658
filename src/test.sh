#!/bin/bash

#test.sh -- batch testing for programs
#usage: ./test.sh <program_path> <tests_dir>

#path to executable
program=$1
#directory containing all .in and .res
test_dir=$2
#levels of verbosity:
# 0: print nothing.
# 1: print diff when wrong.
# 2: print output/expected output when wrong in addition to level 1 actions.
# 3: always print output/expected output in addition to levels 1, 2 actions.
verbose=$3

#setting/validating input
[ -z "$program" ] && { echo "error: must specify executable path"; exit 1; }
[ -z "$test_dir" ] && test_dir="./tests"
[ -z "$verbose" ] && verbose=1
([ ! -x "$program" ] || [ ! -f "$program" ]) && \
    { echo "error: '$program' not a executable file!"; exit 1; }
[ ! -d "$test_dir" ] && { echo "error: '$test_dir' not a directory!"; exit 1; }

#checking if required commands exist
for c in "diff" "basename" "find" "sort" "cut"; do
    command -v $c >/dev/null 2>&1 ||\
        { echo >&2 "error: '$c' program required but not installed."; exit 1; }
done

#counter for passed tests
yes=0
#counter for wrong tests
no=0
#command prefix.
cmd=""
#cmd+="valgrind --leak-check=yes "

#temporary file to store outputs
out="cdmEja22Jjaljb39duh9k.tmp"
#out=$(mktemp)
[ -f "$out" ] && { echo "error: '$out' exists, delete it to proceed."; exit 1; }

#determining input/output extensions
for in_pattern in ".in" ".in.txt" ""; do
    [[ ! -z $(find "$test_dir" -name "*$in_pattern") ]] && break;
done
[ -z "$in_pattern" ] && { echo "error: could not find input files"; exit 1; }
for out_pattern in ".res" ".res.txt" ".out" ".out.txt" ""; do
    [[ ! -z $(find "$test_dir" -name "*$out_pattern") ]] && break;
done
[ -z "$out_pattern" ] && { echo "error: could not find output files"; exit 1; }

echo
echo "--- begin of tests ---"
echo

#iterating over tests
for test_case in $(find "$test_dir" -name "*$in_pattern" | sort); do
	printf -- "--- in $test_case --- "

    #running command for input
	{ $cmd./"$program" < $test_case; } &> $out

    #getting difference between output and expected output
    expected="$test_dir"/$(basename "$test_case" | cut -f1 -d. )"$out_pattern"
    df=$(diff "$out" "$expected")

    #updating list of correct/incorrect test cases
	if [[ -z "$df" ]]; then
        yes=$((yes+1))
        corrects+=$(basename "$test_case")" "
        echo "OK"
    else
        no=$((no+1))
        incorrects+=$(basename "$test_case")" "
        echo "INCORRECT"
    fi

    #printing info if required
    if [[ $verbose -ge 3 ]]; then
        echo "<<<< program output: <<<<"
        cat -- "$out"
        echo ">>>> expected output: >>>>"
        cat "$expected"
    fi
    if [[ ! -z "$df" ]]; then
        if [[ $verbose -eq 2 ]]; then
            echo "<<<< program output: <<<<"
            cat -- "$out"
            echo ">>>> expected output: >>>>"
            cat "$expected"
        fi
        if [[ $verbose -ge 1 ]]; then
            echo "<<<< diff result: >>>>"
            echo "$df"
        fi
    fi
done

#removing temp file
rm -f "$out"

echo
echo "--- end of tests ---"
echo

#summing up results
echo "correct results ($yes):"
echo -e "\t $corrects"
echo "incorrect results ($no):"
echo -e "\t $incorrects"
