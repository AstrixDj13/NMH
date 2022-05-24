% Course: NMH-Lab  
% Semester: SoSe2021
%
% Assignment 3 
%
% This is a Matlab code to solve 1D steady advection-diffusion equation
% discritised by finite-volume schemes 
%
% differential form of advection-diffusion eqn.:
%
%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2
%
% You are asked to fill in the missing parts to complete the implementation.
% Missing parts are marked by ???
%
% Author: Yoshiyuki Sakai
% Email: yoshiyuki.sakai@tum.de


% Clear all variables and plots.
format long;
clear;
clc;
close all

% Set User input for variables and for how many different grid
% resolutions the program should run
prompt = ["Give a value for advection velocity:", "Give a value for diffusivity:", "How many different grid resolutions?:"];
dlgtitle = 'Set variables and cell number';
dims = [1 45];
definput = {'1.0','1.0','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput,'on');

% Let's the user insert the number of grid cells 
for i= 1: str2double(answer{3}(1,:))
    prompt = ["Insert number of grid cells"];
    dlgtitle = 'Set number of cells for different resolutions';
    dims = [1 45];
    grid{i} = inputdlg(prompt,dlgtitle,dims);
end

% Set convection velocity
U0 = str2double(answer{1}(1,:));

% Set diffusivity
Gamma = str2double(answer{2}(1,:));

% Set up grid cells
xend = 2.0 * pi;

% Define figure for plot and set in invisible
f1 = figure("Name", "Comparison of numerical and exact solution with finite volume","Visible","off");

% Initialise array for plot legend before for loop
legend_array = {};

%% Calculation of Numerical and exact solution

for j=1:length(grid)

% Get cell number from user input
cells = str2double(grid{j}(1,:)); 

% Define spacing dx
dx = xend/cells;

%Calculate Peclet number
Peclet= U0*dx/Gamma

% Array of grid cell centre locations:
x = dx/2:dx:xend-dx/2;

% Initialization of cell-averaged field
phi = zeros(cells,1);

% Initialization of matrix A
A = zeros(cells,cells);

% Initialization of vector b
b = zeros(cells,1);

% Boundary cell face values
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid cells
% NB: boundary cells are excluded
for i = 2 : cells-1

     a_w = (U0/2 + Gamma/dx);
     a_p = -2*Gamma/dx;
     a_e = ((-U0/2) + Gamma/dx);
     
%     assign values to LHS matrix A
     A(i,i) = a_p;
     A(i,i-1) = a_w;
     A(i,i+1) = a_e;

%     assign values to RHS vector b
     b(i) = 0;
end

% Boundary conditions (Dirichlet at boundary cell faces)

% at i = 1
 A(1,1) = ((U0/2) - (3*Gamma/dx));
 A(1,2) = ((-U0/2) + (Gamma/dx));
 b(1) = -(U0 + (2*Gamma/dx))*phi_0;

% at i = cells
 A(cells,cells) = ((U0/2) - (3*Gamma/dx));
 A(cells,cells-1) = ((U0/2) + (Gamma/dx));
 b(cells) = (U0 - (2*Gamma/dx))*phi_end;

% Solution of the linear system and save it for each resolution
phi = A\b;
phi_grid{j,1}=phi;
phi_grid{j,2}=x;

% Compute analytical solution
phi_analytic = ((exp((U0*x)/Gamma))-1)/((exp((2*pi*U0)/Gamma))-1);

% Compute relative error at x=pi
nn = ceil(cells / 2);
er(j,2) = abs((phi_analytic(nn) - phi(nn).')/(phi_analytic(nn)));
er(j,1) = dx;

% Compute mean error
err_mean(j,2) = sqrt(mean((phi_analytic - phi.').^2))/mean(phi_analytic);
err_mean(j,1) = dx;

% Plot the numerical and analytical solutions
set(0, 'CurrentFigure',f1);
hold all; 
plot(phi_grid{j,2},phi_grid{j,1},'Linewidth',2);
legend_array{j} = "Number of points = " + grid{j};

end
%% PLOTTING 
   
% Plot the exact solution for the finest availaible grid resolution and
% show figure
  plot(x,phi_analytic,'-go','Linewidth',2);
  legend_array{end+1} = "Exact";
  legend(legend_array,'location','Southeast');
  set(gca,'FontSize',14); 
  f1.WindowState = 'maximized';
  set(f1,'Visible','on');
  xlabel('x');
  ylabel('Phi');
  title("Finite Volume");

% Plot relative error
  f2 = figure("Name","Error plot log scale");
  f2.WindowState = 'maximized';
  loglog(er(:,1),er(:,2),'-bo',x,x.^2,x,x,'LineWidth',2);
  title('Relative Error plot log scale');
  xlabel('Grid spacing dx');
  ylabel('Relative error');
  set(gca,'FontSize',14); 
  legend('Relative Error','Quadratic function','Linear','Location','Southeast')

  % Plot mean error
  f3 = figure("Name","Error plot log scale");
  f3.WindowState = 'maximized';
  loglog(err_mean(:,1),err_mean(:,2),'-bo',x,x.^2,x,x,'LineWidth',2);
  title('Mean Error plot log scale');
  xlabel('Grid spacing dx');
  ylabel('Mean error');
  set(gca,'FontSize',14); 
  legend('Mean Error','Quadratic function','Linear','Location','Southeast')