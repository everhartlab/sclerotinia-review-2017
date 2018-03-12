#!/bin/bash

CITELINE=$(grep -n ^.section..References review.tex | awk -F: '{print $1}')
head -n $(( CITELINE - 1 )) review.tex > tmp.tex
cat tmp.tex formatted-citations.tex > review.tex
rm tmp.tex
