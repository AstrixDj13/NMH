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
clear;
hold off;

% Set convection velocity
U0 = 1.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 100; 
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;
tend   = dt * tsteps;
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

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

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

  % Analytical solution
   t=i*dt;
  for j = 1 : points
    phi_a(j) = sin(x(j)-U0*t);
    
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.01);

end

