#!/bin/bash
# Use this script to remove all run.sh and setup.sh generated files and folders

# remove log files
rm *_log.txt
# remove folders 1, 3-7
rm -r [1,3-7]*
