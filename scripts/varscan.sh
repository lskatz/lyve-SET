#!/bin/bash

# script to call varscan.jar
# adapted from Seth Sims <xzy3@cdc.gov> September 26, 2014
# Lee Katz <lkatz@cdc.gov>

JAR=$(dirname $0)"/../lib/varscan.v2.3.7.jar"
if [ ! -e "$JAR" ]; then
  echo "ERROR: jar file does not exist at $JAR"
  echo "  Please edit $0 to reflect where the jar file is."
  exit 1;
fi

if [ "$1" == "" ]; then
  script=$(basename $0)
  echo "This is the Lyve-SET method of calling varscan. Usage help below is from varscan itself."
  echo "Run 'java -jar $JAR' for more information"
  echo
  echo "USAGE: $script [COMMAND] [OPTIONS]"
  echo
  java -jar "$JAR" 2>&1 | grep -v 'java -jar' | grep .
  echo
  exit 1
fi

java -jar "$JAR" "$@"
