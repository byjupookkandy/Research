#!/bin/bash
cat sst-hmix-clim.ps | ps2eps -B -s a4 > test.eps
convert +antialias -density 720 test.eps -background  White -flatten sst-hmix-clim.png
#cat ssta-spectrum-chunks-5000yrs.eps | ps2eps -B -s a4 > test.eps
convert +antialias -density 720 ssta-spectrum-chunks-5000yrs.eps -background  White -flatten ssta-spectrum-chunks-5000yrs.png
convert +antialias -density 720 ssta-spectrum-chunks-2K-1K-500yr.eps -background  White -flatten ssta-spectrum-chunks-2K-1K-500yr.png
