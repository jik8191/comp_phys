#!/bin/bash

# need to run bibtex to get references right

pdflatex writeup.tex
bibtex writeup
pdflatex writeup.tex
pdflatex writeup.tex
