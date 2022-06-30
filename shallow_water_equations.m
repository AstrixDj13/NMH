function [ dhdt , dhudt, dhvdt ] = shallow_water_equations( ...
    constants, grid, t, h, hu, hv, kst, zb, dhdt, dhudt, dhvdt)
%SHALLOW_WATER_EQUATIONS Evaluate the time derivatives in the shallow water
%equations

%% References
% We use a Finite-Volume scheme with the approximate Riemann solver from
%
% P. Roe. Approximate Riemann Solvers, Parameter vectors, and Difference
% Schemes, Journal of Computational Physics, 43(2), pp. 357-372, 1981.
%
% P. Glaister. A weak formulation of Roe's approximate Riemann solver
% applied to the St. Venant equations. J. Comput. Phys., 116:189-191, 1995.
%
% www-cfd.ifh.uni-karlsruhe.de/uhlmann/reports_comp/shallow/report_12.html

% Reset the time derivative fields
dhdt(:,:) = NaN;
dhudt(:,:) = NaN;
dhvdt(:,:) = NaN;

% Cell area
A = grid.dx * grid.dy;

%% Preallocate for speed
% Create vectors
alpha_e = zeros(3,1);
alpha_w = zeros(3,1);
alpha_n = zeros(3,1);
alpha_s = zeros(3,1);

beta_e = zeros(3,1);
beta_w = zeros(3,1);
beta_n = zeros(3,1);
beta_s = zeros(3,1);

lambda_e = zeros(3,1);
lambda_w = zeros(3,1);
lambda_n = zeros(3,1);
lambda_s = zeros(3,1);

% Create matrices
V_e = zeros(3,3);
V_w = zeros(3,3);
V_n = zeros(3,3);
V_s = zeros(3,3);

% Loop over interior grid cells
for j = grid.NGHOST+1:length(grid.y)-grid.NGHOST    
    % The i-loop is the innermost loop because in MATLAB arrays are stored
    % in column-major order
    for i = grid.NGHOST+1:length(grid.x)-grid.NGHOST
        %% Compute values at cell faces
        % As explained by Roe (1981), for the shallow water equation one
        % can define a set of variables
        % w_1 = sqrt(h)
        % w_2 = sqrt(h)*u
        % w_3 = sqrt(h)*v
        % such that all fluxes are quadratic in those variables. From those
        % variables, one can derive how to define the means at the cell
        % faces to evaluate the flux Jacobians.
        
        % East face
        h_e = 0.5 * ( h(i,j) + h(i+1,j) );
        u_e = ( hu(i,j) / sqrt( h(i,j) ) + hu(i+1,j) / sqrt( h(i+1,j) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i+1,j) ) );
        v_e = ( hv(i,j) / sqrt( h(i,j) ) + hv(i+1,j) / sqrt( h(i+1,j) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i+1,j) ) );
        
        % West face
        h_w = 0.5 * ( h(i,j) + h(i-1,j) );
        u_w = ( hu(i,j) / sqrt( h(i,j) ) + hu(i-1,j) / sqrt( h(i-1,j) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i-1,j) ) );
        v_w = ( hv(i,j) / sqrt( h(i,j) ) + hv(i-1,j) / sqrt( h(i-1,j) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i-1,j) ) );
        
        % North face
        h_n = 0.5 * ( h(i,j) + h(i,j+1) );
        u_n = ( hu(i,j) / sqrt( h(i,j) ) + hu(i,j+1) / sqrt( h(i,j+1) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i,j+1) ) );
        v_n = ( hv(i,j) / sqrt( h(i,j) ) + hv(i,j+1) / sqrt( h(i,j+1) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i,j+1) ) );
        
        % South face
        h_s = 0.5 * ( h(i,j) + h(i,j-1) );
        u_s = ( hu(i,j) / sqrt( h(i,j) ) + hu(i,j-1) / sqrt( h(i,j-1) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i,j-1) ) );
        v_s = ( hv(i,j) / sqrt( h(i,j) ) + hv(i,j-1) / sqrt( h(i,j-1) ) ) / ( sqrt( h(i,j) ) + sqrt( h(i,j-1) ) );
        
        %% Compute eigensystem
        % Compute wave speeds
        c_e = sqrt( constants.g * h_e );
        c_w = sqrt( constants.g * h_w );
        c_n = sqrt( constants.g * h_n );
        c_s = sqrt( constants.g * h_s );
        
        % Right eigenvectors
        V_e(1,1) = 0.0;
        V_e(2,1) = 0.0;
        V_e(3,1) = 1.0;
        V_e(1,2) = 1.0;
        V_e(2,2) = u_e + c_e;
        V_e(3,2) = v_e;
        V_e(1,3) = 1.0;
        V_e(2,3) = u_e - c_e;
        V_e(3,3) = v_e;
        
        V_w(1,1) = 0.0;
        V_w(2,1) = 0.0;
        V_w(3,1) = 1.0;
        V_w(1,2) = 1.0;
        V_w(2,2) = u_w + c_w;
        V_w(3,2) = v_w;
        V_w(1,3) = 1.0;
        V_w(2,3) = u_w - c_w;
        V_w(3,3) = v_w;
        
        V_n(1,1) = 0.0;
        V_n(2,1) = 1.0;
        V_n(3,1) = 0.0;
        V_n(1,2) = 1.0;
        V_n(2,2) = u_n;
        V_n(3,2) = v_n + c_n;
        V_n(1,3) = 1.0;
        V_n(2,3) = u_n;
        V_n(3,3) = v_n - c_n;
        
        V_s(1,1) = 0.0;
        V_s(2,1) = 1.0;
        V_s(3,1) = 0.0;
        V_s(1,2) = 1.0;
        V_s(2,2) = u_s;
        V_s(3,2) = v_s + c_s;
        V_s(1,3) = 1.0;
        V_s(2,3) = u_s;
        V_s(3,3) = v_s - c_s;
        
        % Eigenvalues (characteristic speeds)
        lambda_e(1) = u_e;
        lambda_e(2) = u_e + c_e;
        lambda_e(3) = u_e - c_e;
        
        lambda_w(1) = u_w;
        lambda_w(2) = u_w + c_w;
        lambda_w(3) = u_w - c_w;
        
        lambda_n(1) = v_n;
        lambda_n(2) = v_n + c_n;
        lambda_n(3) = v_n - c_n;
        
        lambda_s(1) = v_n;
        lambda_s(2) = v_n + c_s;
        lambda_s(3) = v_n - c_s;
        
        %% Volume, convective and pressure fluxes
        % Compute difference vectors
        diff_e = [h(i+1,j)-h(i,j); hu(i+1,j)-hu(i,j); hv(i+1,j)-hv(i,j)];
        diff_w = [h(i,j)-h(i-1,j); hu(i,j)-hu(i-1,j); hv(i,j)-hv(i-1,j)];
        diff_n = [h(i,j+1)-h(i,j); hu(i,j+1)-hu(i,j); hv(i,j+1)-hv(i,j)];
        diff_s = [h(i,j)-h(i,j-1); hu(i,j)-hu(i,j-1); hv(i,j)-hv(i,j-1)];
        
        % Characteristic decomposition of flux differences
        alpha_e(1) = -v_e * diff_e(1) + diff_e(3);
        alpha_e(2) = ( ( c_e - u_e) * diff_e(1) + diff_e(2) ) / ( 2*c_e );
        alpha_e(3) = ( ( c_e + u_e) * diff_e(1) - diff_e(2) ) / ( 2*c_e ); 
        
        alpha_w(1) = -v_w * diff_w(1) + diff_w(3);
        alpha_w(2) = ( ( c_w - u_w) * diff_w(1) + diff_w(2) ) / ( 2*c_w );
        alpha_w(3) = ( ( c_w + u_w) * diff_w(1) - diff_w(2) ) / ( 2*c_w );
        
        alpha_n(1) = -u_n * diff_n(1) + diff_n(2);
        alpha_n(2) = ( ( c_n - v_n) * diff_n(1) + diff_n(3) ) / ( 2*c_n );
        alpha_n(3) = ( ( c_n + u_n) * diff_n(1) - diff_n(3) ) / ( 2*c_n ); 

        alpha_s(1) = -u_s * diff_s(1) + diff_s(2);
        alpha_s(2) = ( ( c_s - v_s) * diff_s(1) + diff_s(3) ) / ( 2*c_s );
        alpha_s(3) = ( ( c_s + u_s) * diff_s(1) - diff_s(3) ) / ( 2*c_s );
 
        % Volume fluxes
        q_e = 0.5 *( hu(i,j) + hu(i+1,j) - V_e(1,:) * ( abs( lambda_e ) .* alpha_e ) ) * grid.dy;
        q_w = 0.5 *( hu(i,j) + hu(i-1,j) - V_w(1,:) * ( abs( lambda_w ) .* alpha_w ) ) * grid.dy;
        q_n = 0.5 *( hv(i,j) + hv(i,j+1) - V_n(1,:) * ( abs( lambda_n ) .* alpha_n ) ) * grid.dx;
        q_s = 0.5 *( hv(i,j) + hv(i,j-1) - V_s(1,:) * ( abs( lambda_s ) .* alpha_s ) ) * grid.dx;
        
        % Momentum fluxes
        qu_e = 0.5 *( hu(i,j)^2 / h(i,j) + 0.5 * constants.g * h(i,j) ^2 ...
            + hu(i+1,j)^2/ h(i+1,j) + 0.5 * constants.g * h(i+1,j)^2 ...
            - V_e(2,:) * ( abs( lambda_e ) .* alpha_e ) ) * grid.dy;
        qu_w =  0.5 *( hu(i,j)^2 / h(i,j) + 0.5 * constants.g * h(i,j) ^2 ...
            + hu(i-1,j)^2/ h(i-1,j) + 0.5 * constants.g * h(i-1,j)^2 ...
            - V_w(2,:) * ( abs( lambda_w ) .* alpha_w ) ) * grid.dy;
        qu_n = 0.5 *( hu(i,j) * hv(i,j) / h(i,j) + hu(i,j+1) * hv(i,j+1) / h(i,j+1) ...
            - V_n(2,:) * ( abs( lambda_n ) .* alpha_n ) ) * grid.dx;
        qu_s = 0.5 *( hu(i,j) * hv(i,j) / h(i,j) + hu(i,j-1) * hv(i,j-1) / h(i,j-1) ...
            - V_s(2,:) * ( abs( lambda_s ) .* alpha_s ) ) * grid.dx;
        
        qv_e = 0.5 *( hu(i,j) * hv(i,j) / h(i,j) + hu(i+1,j) * hv(i+1,j) / h(i+1,j) ...
            - V_e(3,:) * ( abs( lambda_e ) .* alpha_e ) ) * grid.dy;
        qv_w = 0.5 *( hu(i,j) * hv(i,j) / h(i,j) + hu(i-1,j) * hv(i-1,j) / h(i-1,j) ...
            - V_w(3,:) * ( abs( lambda_w ) .* alpha_w ) ) * grid.dy;
        qv_n = 0.5 *( hv(i,j)^2 / h(i,j) + 0.5 * constants.g * h(i,j) ^2 ...
            + hv(i,j+1)^2/ h(i,j+1) + 0.5 * constants.g * h(i,j+1)^2 ...
            - V_n(3,:) * ( abs( lambda_n ) .* alpha_n ) ) * grid.dx;
        qv_s = 0.5 *( hv(i,j)^2 / h(i,j) + 0.5 * constants.g * h(i,j) ^2 ...
            + hv(i,j-1)^2/ h(i,j-1) + 0.5 * constants.g * h(i,j-1)^2 ...
            - V_s(3,:) * ( abs( lambda_s ) .* alpha_s ) ) * grid.dx;
        
        
        %% Bottom slope
        % Coefficients for characteristic decomposition of source term     
        beta_e(1) = 0.0;
        beta_e(2) = 0.25 * ( sign( lambda_e(2) ) - 1 ) * c_e *  ( zb(i+1,j) - zb(i,j) );
        beta_e(3) = 0.25 * ( 1 - sign( lambda_e(3) ) ) * c_e *  ( zb(i+1,j) - zb(i,j) );
        
        beta_w(1) = 0.0;
        beta_w(2) = 0.25 * ( sign( lambda_w(3) ) - 1 ) * c_w *  ( zb(i,j) - zb(i-1,j) );
        beta_w(3) = 0.25 * ( 1 - sign( lambda_w(2) ) ) * c_w *  ( zb(i,j) - zb(i-1,j) );
        
        beta_n(1) = 0.0;
        beta_n(2) = 0.25 * ( sign( lambda_n(2) ) - 1 ) * c_n *  ( zb(i,j+1) - zb(i,j) );
        beta_n(3) = 0.25 * ( 1 - sign( lambda_n(3) ) ) * c_n *  ( zb(i,j+1) - zb(i,j) );
        
        beta_s(1) = 0.0;
        beta_s(2) = 0.25 * ( sign( lambda_s(3) ) - 1 ) * c_s *  ( zb(i,j) - zb(i,j-1) );
        beta_s(3) = 0.25 * ( 1 - sign( lambda_s(2) ) ) * c_s *  ( zb(i,j) - zb(i,j-1) );

        % Compute source term in eigenvector basis
        S_bottom = V_e * beta_e * grid.dy ...
            + V_w * beta_w * grid.dy ...
            + V_n * beta_n * grid.dx ... 
            + V_s * beta_s * grid.dx;

        %S_bottom = [0;-h(i,j)*constants.g*(zb(i+1,j)-zb(i-1))/(2*grid.dx)*A;-h(i,j)*constants.g*(zb(i+1,j)-zb(i-1))/(2*grid.dx)*A]
        
        %% Bottom shear stress
        q = sqrt( hu(i,j)^2 + hv(i,j)^2 );
        
        % x-direction
        s_friction_x = - A * constants.g * q * hu(i,j) / ( kst(i-grid.NGHOST,j-grid.NGHOST)^2 * h(i,j)^(4/3) );
        
        % y-direction
        s_friction_y = - A * constants.g * q * hv(i,j) / ( kst(i-grid.NGHOST,j-grid.NGHOST)^2 * h(i,j)^(4/3) );
        
        %% Summing up
        % Continuity equation
        dhdt(i,j) = ( -q_e + q_w - q_n + q_s + S_bottom(1) ) / A ;
        
        % u-momentum equation
        dhudt(i,j) = ( - qu_e + qu_w - qu_n + qu_s ...
            + S_bottom(2) + s_friction_x ) / A;
        
        % v-momentum equation
        dhvdt(i,j) = ( - qv_e + qv_w - qv_n + qv_s ...
            + S_bottom(3) + s_friction_y ) / A;
        
    end
end


end

