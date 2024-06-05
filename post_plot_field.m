clc; clear
addpath('postprocessFunctions')
folderName   = 'resFolder';
scaleFactor  = 0;
varName  = 'dmg';
step     = 500;
gridSize = 8;
plotGivenFieldVariables(folderName, scaleFactor, varName, step, gridSize)