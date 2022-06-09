% Course: CFD Lab
% TU-Muenchen, SS 2018
%
% This is the Matlab script for the unsteady 1D convection diffusion equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% F. Mintgen, f.mintgen@bv.tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi             d^2 phi
%  ------- =  -U0 * ------- + Gamma * ---------
%    dt               dx                dx^2 
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
% Crank Nicolson scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear
close all
clc;
hold off;

% Set convection velocity
U0 = 1.0;

% Set diffusion coefficient
Gamma = 1.0;

% Discrete spacing in space
xend   = 2.0 * pi;
points = 40; 
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
dt     = 0.1;
tend   = dt * tsteps;

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

phi_time(1,:) = phi;

% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

% Compute coefficients of matrix A
  a_w = 0.5 * dt *(-U0/2./dx - Gamma/dx^2.);
  a_p = 1. + dt * Gamma / dx^2.;
  a_e = 0.5 * dt *(U0/2./dx - Gamma/dx^2.);

% Crank Nicolson:
%----------------
%
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
    A(1,points-1) = a_w;
    A(1,1) = a_p;
    A(1,2) = a_e; 
    b(1) = phi(points-1) * (-a_w) + phi(1) * (1. - dt * Gamma /dx^2.) + phi(2) * (-a_e);

  % Loop over grid points in space:
  for j = 2 : points - 1

     A(j,j-1) = a_w;
     A(j,j) = a_p; 
     A(j,j+1) = a_e;
     b(j) = phi(j-1) * (-a_w) + phi(j) * (1. - dt * Gamma /dx^2.) + phi(j+1) * (-a_e);

  end

  % Periodic boundary conditions at x=2*pi:
    A(points,2) = a_e;
    A(points,points-1) = a_w;
    A(points,points) = a_p;
    b(points) = b(1);

  % Solve the linear system of equations
  phi = A\b;

  phi_time(i+1,:)=phi;

  % Analytical solution
  t = i*dt;
  for j = 1 : points
    phi_a(j) = exp(-Gamma*t)*sin(x(j)-U0*t);
    
  end

  % Plot transported wave for each timestep
  fg1 = figure(1);
  plot(x, phi, 'r', x, phi_a, 'g');
  hold off;
  pause(0.001);

end

pause(2);
close(fg1);

%% Saving plots

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
            disp("You entered a wrong time step.")
            break;
        else
            phi_analytical = exp(-Gamma*time(i)*dt)*sin(x-U0*dt*time(i));
            fh=figure();
            fh.WindowState = "maximized";
            plot(x,phi_time(time(i)+1,:),'r',x,phi_analytical,'g','LineWidth',2);
            set(gca, 'FontSize',15);
            xlabel('x');
            ylabel('\phi (x,t)');
            legend('Numerical solution','Analytical solution','Location','Southeast');
            title_str = "Crank-Nicolson, Time = " + time(i)*dt + " ,dt = " + dt + " ,Points = " + points;
            title(title_str);
            save_str = "Nicolson_Timestep_" + time(i);
            print('-djpeg',save_str,'-r250');
        end

    end
end