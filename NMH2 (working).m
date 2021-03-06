% Course: CFD-Lab
% SS 2020

% This is the Matlab script for the convection-diffusion equation 

% You must fill in the missing parts by yourself!
% Missing parts are marked by ???

% Tianshi Sun
% tianshi.sun@tum.de


% solution of

%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2

% u=0
% 0 <= x <= 2pi

% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
U0 = -1.0;

% Set diffusivity
Gamma = 1.0;

% Discrete spacing in space!
xend   = 2.0*pi;
points = 501; 
dx     = xend/(points-1);
% Grid with x locations:
x = (0.0 : dx : xend).';

% Initialization of field
phi = zeros(points,1);

% Initialization of matrix A
A = zeros(points,points);

% Initialization of vector b
b = zeros(points,1);

% Boundary condition
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid points in space
% note that boundary points are excluded
for i = 2 : points-1

     a_w = (U0/dx)+Gamma/(dx*dx);
     a_p = -(U0/dx) -(2*Gamma)/(dx*dx);
     a_e = Gamma/(dx*dx);
     
%     assign values to matrix A
      A(i,i-1) = a_w;
      A(i,i) = a_p;
      A(i,i+1) = a_e;

%     assign values to vector b


end

% Boundary conditions

% at i = 1
 A(1,1) = 1;
 b(1) = phi_0;

% at i = points
 A(points,points) = 1;
 b(points) = phi_end;

% Solution of the linear system
phi = A\b;

%Analytical solution

phi_analytic = (exp((U0*x)/Gamma)-1)/(exp((2*pi*U0)/Gamma)-1);

%error
pointIsNotOnPi = mod(points,2);

if pointIsNotOnPi == false
    disp("Please choose an odd grid resolution so that points/2 = π")
end

 % Plot the solution
 figure();
 plot(x,phi,'r', x,phi_analytic, '-go');
 title("Number of points = "+ points);

nn = (points + 1) / 2;
er = abs(phi_analytic(nn) - (phi(nn).')/phi_analytic(nn));

% Plot error not fully implemented
figure();
loglog(er,dx,'-bo')