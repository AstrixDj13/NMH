% Course: CFD Lab
% TU Muenchen, Summer term 2020
%
% This is the Matlab script for the unsteady 1D convection equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% Author: Tianshi Sun, tianshi.sun@tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi
%  ------- =  -U0 * -------
%    dt               dx
%
% 0 <= x <= 2pi
%
% periodic boundary condition
% phi(0) = phi(2pi)
%
% initial condition
% t = 0  ==>  phi = sin(x)
%
% Central difference scheme (CDS) for spatial discretization
% Implicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear all;
close all;
clc;
hold off;

%% Setting up the simulation

% Set convection velocity
U0 = 10.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 40; 
dx     = xend / ( points - 1 );

% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.01;
tend   = dt * tsteps;

% CFL Number
CFL = U0*dt/dx

% Initialise coefficient matrix A, constant vector b
% and solution vector phi
A   = zeros(points,points);
b   = zeros(points,1);
phi = zeros(points,1);

% Initialise the solution (initial condition)
% Loop over grid points in space:
for j = 1 : points
   phi(j) = sin(x(j));
end

% Store the solution for every timestep in phi_time
phi_time(1,:)=phi;

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

%% Calculation using Implicit Euler

% Compute coefficients of matrix A
 a_w = -U0*dt/(2*dx);
 a_p = 1;
 a_e = U0*dt/(2*dx);

% Implicit Euler:
%----------------
%
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
   A(1,1) = a_p;
   A(1,2) = a_e;
   A(1,points-1) = a_w;
   b(1) = phi(1);

  % Loop over grid points in space:
  for j = 2 : points - 1

     A(j,j-1) = a_w;
     A(j,j) = a_p;
     A(j,j+1) = a_e;
     b(j) = phi(j);

  end

  % Periodic boundary conditions at x=2*pi:
   A(points,2) = a_e;
   A(points,points-1) = a_w;
   A(points,points) = a_p;
   b(points) = phi(points);

  % Solve the linear system of equations
  phi = A\b;

  % Add current solution to phi_time
  phi_time(i+1,:)=phi;

  % Analytical solution
   t=i*dt;
  for j = 1 : points
    phi_a(j) = sin(x(j)-U0*t);
    
  end

  % Plot transported wave for each timestep
  fg1=figure(1);
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.001);

end

pause(2);
close(fg1);

%% Saving Plots

prompt = ['Type comma separated the number of time steps for which to save the figure. Range [0-' num2str(tsteps) ']'];
dlgtitle = 'Save';
dims = [1 50];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
time= str2double(split(answer{1},','));

if answer ==""
    disp("Not plots will be saved")
else

    for i=1:length(time)
        if time(i) >tsteps || time(i)<0
            disp("You entered a wrong time step")
            break;
        else
            phi_analytical = sin(x-U0*dt*time(i));
            fh =figure();
            fh.WindowState = "maximized";
            plot(x,phi_time(time(i)+1,:),'r',x,phi_analytical,'g','LineWidth',2);
            set(gca, 'FontSize',15);
            xlabel('x');
            ylabel('\phi (x,t)');
            legend('Numerical solution','Analytical solution','Location','Southeast');
            title_str = "Implicit Euler, Time = " + time(i)*dt + " ,dt = " + dt + " ,Points = " + points;
            title(title_str);
            save_str = "Implicit_Timestep_" + time(i);
            print('-djpeg',save_str,'-r250');
        end

    end
end
