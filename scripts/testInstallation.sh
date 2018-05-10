#!/bin/bash

set -e

set_test.pl --numcpus 2 lambda lambda -- --noqsub

echo "Everything passed!"

