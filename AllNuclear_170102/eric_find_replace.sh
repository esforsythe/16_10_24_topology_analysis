#!/bin/bash

#replace.txt is a file with several find and replace phrases
#example: 's/find/replace/g'

while read p; do
  #echo $p
  sed -ibackup $p file.txt
done <replace.txt