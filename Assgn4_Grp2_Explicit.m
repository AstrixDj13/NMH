%Course: NMH Lab
% TU Muenchen, Summer term 2021
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
% Explicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear all;
close all;
clc
hold off;

% Set convection velocity
U0 = 1.0;

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

CFL = U0*dt/dx

% Initialise the solution (initial condition)
% Loop over grid points in space:
for j = 1 : points
   phi(j) = sin(x(j));
end
phi_time(1,:)=phi;
% Check initial field:
plot(x, phi, 'r');
hold on;
pause(3);

% Explicit Euler:
%----------------
%
% phinew is phi at new time level
% phinew must be written back to phi for new timestep
Aw = U0*dt/(2*dx);
Ap =1;
Ae = -U0*dt/(2*dx);
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
  phinew(1) = Aw*phi(points-1) + Ap*phi(1)+Ae*phi(2);

  % Loop over grid points in space:
  for j = 2 : points - 1

    phinew(j) = Aw*phi(j-1)+Ap*phi(j)+Ae*phi(j+1);

  end

  % Periodic boundary conditions at x=2*pi:
  phinew(points) = Aw*phi(points-1)+Ap*phi(points)+Ae*phi(2);

  % Write new field back to old field:
  phi = phinew;
  phi_time(i+1,:)=phinew;
  % Analytical solution
  
  for j = 1 : points
    t = i*dt;
    phi_a(j) = sin(x(j)-U0*t);
    
  end

  % Plot transported wave for each timestep
  fg1= figure(1);
  plot(x, phi, 'r', x, phi_a, 'g');
  %ylim([-2,2])
  hold off;
  pause(0.001);
  %fprintf('%d\n',i)
end
close(fg1);
%% Saving plots

prompt = "Type the number of time steps to save the figure [0-1000]? Type it comma separated";
dlgtitle = 'Save';
dims = [1 50];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
time= str2double(split(answer{1},','));

if answer ==""
    disp("Not plots will be saved")
else

    for i=1:length(time)
        if time(i) >1000 || time(i)<0
            disp("You entered a wrong time step.")
            break;
        else
            phi_analytical = sin(x-U0*dt*time(i));
            figure()
            plot(x,phi_time(time(i)+1,:),'r',x,phi_analytical,'g');
            legend('Numerical solution','Analytical solution','Location','Southeast');
            title_str = "Explicit Euler, Time = " + time(i)*dt + " ,dt = " + dt + " ,Points = " + points;
            title(title_str);
            save_str = "Explicit Euler at Time " + time(i)*dt;
            print('-djpeg',save_str,'-r250');
        end

    end
end
