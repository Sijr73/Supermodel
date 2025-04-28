#!/bin/bash

rsync -avzru --stats --copy-links -e ssh hpc:~/envirDist/fbadata/* fbadata/
rsync -avzru --stats --copy-links -e ssh hpc:~/envirDist/envirdata/* envirdata/
rsync -avzru --stats --copy-links -e ssh hpc:~/envirDist/crossdata/* crossdata/

rsync -avzru --stats --copy-links -e ssh hpc:~/envirDist/envirDistMILP.iAF1260.Rdata .






