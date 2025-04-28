#!/bin/sh


seq 76 | xargs -n 1 -P 10 -I@ sh -c "Rscript fbaPreCalc.R @ &>logfiles/fbaPreCalc.@.out"
