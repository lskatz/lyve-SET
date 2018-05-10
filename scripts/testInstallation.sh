#!/bin/bash

set -e

set_test.pl --numcpus 2 lambda lambda -- --noqsub --nodiagnose

echo "Everything passed!"

