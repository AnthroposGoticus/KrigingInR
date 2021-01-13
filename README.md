# KrigingInR

Script for Linux(tested)/Windows(not Tested)
For graphical output Rstudio is recommended

Functions.R contains most of the main functions for the other scripts

CalcLoop.R - Main script calculates Kriging over multiple polygons within a research area (e.g. a disturbed geological layer)

  to work correctly you have to upload your...
  
  ...samplepoints to "$PWD/Geodaten/Shapes/KrigingSamplepoints"
  
  ...overall modelarea to "$PWD/Geodaten/Shapes/Model"
  
  ...(if you have a single shape containing the subpolygons) "$PWD/Geodaten/Shapes/Model"
  
  ...(if you have a file for each subpolygon) "$PWD/Geodaten/Shapes/Model/Polygone"
  
