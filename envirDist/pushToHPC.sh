#!/bin/bash

rsync -avz --stats --copy-links --exclude-from=excludeFile.txt -e ssh * hpc:~/envirDist/
