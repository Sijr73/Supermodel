#!/bin/bash



for file in plots/Fig{,S}?.pdf; do
	echo $file
	convert -density 300 -border 1x1 -bordercolor white -fuzz 0% -trim $file jpg/`basename $file .pdf`.jpg
done
