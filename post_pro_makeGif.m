clc; clear
addpath('postprocessFunctions')
gridSize = 10;
resFileDir = 'resFolder';
gifFileName = 'omega3';
scaleFactor = 5;
varName = 'm31';
startStep = 10;
step    = 10;
endStep = 490;
makeGif(resFileDir, gifFileName, scaleFactor, varName, startStep, step, endStep, gridSize)