#!/bin/bash -l
grep -v 'Orthogroup' $1 | tr ',' '\t' > $2
cat proteome/*.faa| grep '>' | sed 's/.*>//g' | sed s/' '/';'/1 > $3
