clear all
clc
close all

%Define vector with number of grid points
grid = [10,25,50,100,1000,10000];

for exp = 1:1:6

    % Number of grid points
    n = grid(exp);
    % Grid spacing
    h = 2.0*pi / (n-1);

    % Calculation of... 
    for i = 1 : n

      % ...coordinates
      x(i) = 0.0 + (i-1) * h;

      % ...analytical values of function and derivatives
      f(i) = sin(2*x(i));
      dfe(i) = 2*cos(2* x(i) );
      dffe(i) = -4*sin(2*x(i));

    end

    for i = 1 : n

      % Periodic boundary conditions... 

      % ...for the first point
      if ( i == 1 )

        fp = f(i+1);
        fm = f(n-1);

      % ...for the last point
      elseif ( i == n )

        fp = f(2);
        fm = f(n-1);

      % If point is not at the boundaries 
      else
        fp = f(i+1);
        fm = f(i-1);
      end

      fi = f(i);

    % Numerical approximation of first derivative...
      dfn_up(i) = ( fi - fm ) / ( 1.0*h );
      dfn_down(i) = ( fp - fi ) / ( 1.0*h );
      dfn_central(i) = ( fp - fm ) / ( 2.0*h );

     % ...and second derivative
      dffn_central(i) = (fp - 2*fi + fm)/(h*h);

     % Linear Interpolation 
      interpol(i) = (fi + fp)/2; 

    % Calculation of the error for first derivative...
      er_up(i) = abs( ( dfe(i) - dfn_up(i) ) / dfe(i) );
      er_down(i) = abs( ( dfe(i) - dfn_down(i) ) / dfe(i) );
      er_central(i) = abs( ( dfe(i) - dfn_central(i) ) / dfe(i) );

     % ...and for second derivative
      er_second(i) = abs( ( dffe(i) - dffn_central(i) ) / dffe(i) );

     % and for linear interpolation
      er_interpol(i) = abs((fi - interpol(i) )/fi);

    end

    % Plotting of analytical solution and numerical approximation of first
    % derivative
    f1=figure(1);
    f1.WindowState = 'maximized';
    subplot(3,2,exp)
    plot(x,dfn_up,x, dfn_down,x, dfn_central, x, dfe )
    legend('Upwind','Downwind','Central','Exact','Location','SouthEast')
    set(gca,'FontSize',14);
    title(n);
    pause(1);
    % save plots
  

    % Plotting of analytical solution and numerical approximation of
    % second derivative
    f2=figure(2);
    f2.WindowState = 'maximized';
    subplot(3,2,exp)
    plot(x,dffn_central, x, dffe )
    legend('Central Second Derivative', 'Exact Second Derivative' ,'Location','SouthEast')
    set(gca,'FontSize',14); 
    title(n);
    pause(1);

     % Plotting of analytical solution and numerical approximation of
     % linear interpolation
    f3=figure(3);
    f3.WindowState = 'maximized';
    subplot(3,2,exp)
    plot(x,interpol, x, f )
    legend('Linear Interpolation', 'Exact' ,'Location','SouthEast')
    set(gca,'FontSize',14); 
    title(n);
    pause(1);
    


    % Storing grid spacing and error in matrix 'error' for...
    % first derivative...
    error_up(exp,1) = h;
    error_up(exp,2) = er_up(n/5);

    error_down(exp,1) = h;
    error_down(exp,2) = er_down(n/5);

    error_central(exp,1) = h;
    error_central(exp,2) = er_central(n/5);

    % ... and for second derivative...
    error_second(exp,1) = h;
    error_second(exp,2) = er_second(n/5);

    % and for linear interpolation
    error_interpol(exp,1) = h;
    error_interpol(exp,2) = er_interpol(n/5);

end

% Plotting error function with linear and quatratic function as comparison

f4= figure(4);
f4.WindowState = 'maximized';
loglog(error_up(:,1),error_up(:,2),'-bo',error_down(:,1),error_down(:,2),'-go',error_central(:,1),error_central(:,2),'-ro',error_second(:,1),error_second(:,2),'-ko',error_interpol(:,1),error_interpol(:,2),'-mo',x,x,x, x.^2);
set(gca,'FontSize',14); 
title('Error plot log scale');
xlabel('Grid spacing h');
ylabel('Relative error');
legend('Error Upwind','Error Downwind','Error Central','Error Central Second Derivative','Error Linear Interpolation','Linear Function', 'Quadratic function','Location','SouthEast')

