%**************************************************************************
% CFD Lab Summer Semester 2020, Assignment 6
%
% This code solves the following 2D scalar advection-diffusion problem 
%
% du/dt = -(u0*du/dx + v0*du/dy) + nu*(d^2u/dx^2 + d^2 u/dy^2) 
%
% with central difference spatial discretisation and explicit Euler 
% time integration.
% Periodic boundary condition is applied both in x- and y-direction
%
% authors: D.Quosdorf & Y.Sakai
% May, 2018
%**************************************************************************


% close figures, clear command window and memory
close all; clc; clear

% read infile for automatic input
infilename = 'infile_2d_condiff_scal.mat';

% build structures 'parm' and 'flow'
[parm, flow] = build_structs;

% fill structs 'parm' and 'flow' with data from infile
[parm, flow] = set_params(parm, flow, infilename);

% initialization of flow field
[parm, flow] = initialize(parm, flow);

% start time integration
[parm, flow] = timeint(parm, flow);

%**************************************************************************

 %plot results (Note: please refer to the previous assigments)
 %???