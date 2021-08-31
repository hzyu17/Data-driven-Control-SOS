% gEDME to approximate the infinitesimal of Koopman Operator


close all
clear all
clc
% ------------------------
% Load sampling data 
% ------------------------
load('sampling_data.mat')

%%
% ------------------------------
% Variables and dimensions 
% ------------------------------
pvar x1 x2 u
x = [x1;x2];
nx = length(x); % variable dimension
m = length(u); % input dimension

% ---------------------------------------
% Setting bx as a known polynomial
% ---------------------------------------

% Choice 1. local controller design using LQR, CLF for bx
% variable definition
x = sym('x',[2,1],'real');
u = sym('u',[0,1],'real');
dx = sym('dx',[2,1],'real');

% model, dxdt = fx + gx*u
fx = [x(2); (1 - x(1)^2) * x(2) - x(1)];
gx = [0; 1];
