% Preamble to all scripts.  Clean slate, and set up paths.
clear; clc;
% close all;

addpath calculations;
addpath calculations/boundsgens;
addpath calculations/evgens;
addpath calculations/nodegens;
addpath calculations/polygen;
addpath calculations/std;

addpath plotting;

addpath tests;

% Set up computation environment so that calculations run the same locally
% as well as in parallel.
preamble_numeric()
