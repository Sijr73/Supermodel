#!/bin/bash

rsync -avzru --stats --copy-links -e ssh hpc:~/energyGeneratingCycleRemoval/results/* results/










