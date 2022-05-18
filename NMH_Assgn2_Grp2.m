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
clc;

% Set User input for variables, scheme and for how many different grid
% resolutions the program should run
prompt = ["Which scheme? Type 'central' or 'upwind':","Give a value for convection diffusivity:", "Give a value for diffusivity:", "How many different grid resolutions?:"];
dlgtitle = 'Set variables, scheme and point number';
dims = [1 45];
definput = {'central','10','1.0','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput,'on');

% Let's the user insert the number of grid points 
for i= 1: str2double(answer{4}(1,:))
    prompt = ["Insert number of grid points"];
    dlgtitle = 'Set number of points for different resolutions';
    dims = [1 45];
    grid{i} = inputdlg(prompt,dlgtitle,dims);
end


% Set scheme
scheme = answer{1}(1,:);
% Check if the scheme input is valid
if (scheme ~= "central") && (scheme ~= "upwind")
    disp("Use central or upwind scheme for first derivative!");
    return;
end

% Set convection velocity
U0 = str2double(answer{2}(1,:));
% Set diffusivity
Gamma = str2double(answer{3}(1,:));

% Discrete spacing in space!
xend   = 2.0*pi;

for j=1:length(grid)

    points = str2double(grid{j}(1,:));

    % Check if it is an uneven number, so we get later the error plot at
    % x=pi
    if mod(points,2) == 0
        disp("Please run again and choose an uneven number of points!");
        return;
    end
    
    dx     = xend/(points-1);
    % Grid with x locations:
    x = 0.0 : dx : xend;
    
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
    if scheme == "upwind"
        for i = 2 : points-1
    
         a_w = (U0/dx)+Gamma/(dx*dx);
         a_p = -(U0/dx) -(2*Gamma)/(dx*dx);
         a_e = Gamma/(dx*dx);
         
    %     assign values to matrix A
          A(i,i-1) = a_w;
          A(i,i) = a_p;
          A(i,i+1) = a_e;
          
        end
    elseif scheme == "central"
    
         for i = 2 : points-1
    
         a_w = (U0/2*dx)-Gamma/(dx*dx);
         a_p = -(2*Gamma)/(dx*dx);
         a_e = -((U0/2*dx)-Gamma/(dx*dx));
         
    %     assign values to matrix A
          A(i,i-1) = a_w;
          A(i,i) = a_p;
          A(i,i+1) = a_e;
    
         end
        
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
    
    phi_analytic = ((exp((U0*x)/Gamma))-1)/((exp((2*pi*U0)/Gamma))-1);
    
    % error
    
    nn = (points + 1) / 2;
    er(j,2) = abs(phi_analytic(nn) - (phi(nn).')/phi_analytic(nn));
    er(j,1) = dx;
    
    
    % Plot the solution
    figure();
    plot(x,phi,'r', x,phi_analytic, '-go');
    title("Number of points = "+ points);

end

% Plot error
  figure();
  loglog(er(:,1),er(:,2),'-bo')
    