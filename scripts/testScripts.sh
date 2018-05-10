#!/bin/bash

# This test assumes that set_test.pl has already been run
# and that it generated project directory lambda.

# Run every shell script in unitTests
testsDir=$(dirname $0)/unittests
echo "Running tests in $testsDir"
FAIL=0;
for i in $testsDir/*.sh; do 
  b=$(basename $i)
  echo -n "Running $b .. "
  bash $i

  # If there is any failure, mark it and move on
  if [ $? -gt 0 ]; then
    echo -n "$b failed"
    FAIL=$(($FAIL + 1))
  else
    echo -n "$b passed"
  fi
  echo
done

# Exit with how many tests failed
if [ "$FAIL" -gt 0 ]; then
  echo "Failed $FAIL tests"
  exit $FAIL
fi

