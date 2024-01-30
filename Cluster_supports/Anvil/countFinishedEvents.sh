#!/usr/bin/env bash

for Folder in `ls `
do
    goodEvents=`cat $Folder/job.o* 2>/dev/null | grep "good" | wc -l`
    badEvents=`cat $Folder/job.o* 2>/dev/null | grep "skip" | wc -l`
    echo "$Folder: $goodEvents/$badEvents"
done
