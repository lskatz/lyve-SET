#!/bin/bash

# This test assumes that set_test.pl has already been run
# and that it generated project directory lambda.

# Run every shell script in unitTests
testsDir=$(dirname $0)/unittests
echo "Running tests in $testsDir"
FAIL=0;
for i in $testsDir/*.sh; do 
  b=$(basename $i)
  echo "Running $b"

  # While running the program, indent any output.
  bash $i 2>&1 | perl -lane 'chomp; print "  $_";'
  # Checking the exit code is bash-specific and so the #! line is important
  exit_code=${PIPESTATUS[0]}

  # If there is any failure, mark it and move on
  if [ "$exit_code" -gt 0 ]; then
    echo "$b failed"
    FAIL=$(($FAIL + 1))
  else
    echo "$b passed"
  fi
done

# Exit with how many tests failed
if [ "$FAIL" -gt 0 ]; then
  echo "Failed $FAIL tests"
  exit $FAIL
fi

