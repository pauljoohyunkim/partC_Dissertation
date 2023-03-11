#!/bin/bash
#
texcount -inc -html -v -sum dissertation.tex > results.html
firefox results.html
