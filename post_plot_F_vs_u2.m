clc; clear; 
addpath('postprocessFunctions')
figure(1)
clf
gridSize = 10;
figure(1)
[u1, u2, fx, fy] = caculateReactionForce('resFolder',10, 10, 500);
plot(u2, -fy,'k-.','linewidth',2);
hold on