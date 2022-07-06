function [ lambda_fix ] = entropy_fix_y(V, lambda, alpha, g, ...
        h_L, hu_L, hv_L, h_R, hu_R, hv_R)
%ENTROPY_FIX_Y Compute a modified set of eigenvalues according to the entropy
%fix by LeVeque. This fix prevents a rarefaction shock when the flow state
%changes from supercritical to subcritical.

%% Notes on the input data
% Characteristic decomposition of the difference between the states
%    [h_R;hu_R;hv_R]-[h_L;hu_L;hv_L] = V * alpha


%% References
% M. Pelanti, L. Quartapelle, and L. Vigevano, A review of entropy fixes as
% applied to Roe’s linearization, Politecnico di Milano, 2000.
% https://perso.ensta-paris.fr/~pelanti/ef_PQV.pdf
%
% R. J. LeVeque. Numerical Methods for Conservation Laws, Lectures in
% Mathematics, ETH Zürich, Birkhäuser, Basel, 1990. Chap. 14, Sec. 2.2

%% Allocate variables
% Result: corrected characteristic speeds
lambda_fix = zeros(3,1);

% Characteristic speeds of intermediate states
lambda_L_int = zeros(3,1);
lambda_R_int = zeros(3,1);


%% Computation
% Loop over eigenvalues and eigenvectors
for k = 1:3
    % Construct intermediate states
    h_L_int  = h_L  + V(1,1:k-1) * alpha(1:k-1);
    hu_L_int = hu_L + V(2,1:k-1) * alpha(1:k-1);
    hv_L_int = hv_L + V(3,1:k-1) * alpha(1:k-1);

    h_R_int  = h_L_int  + V(1,k) * alpha(k);
    hu_R_int = hu_L_int + V(2,k) * alpha(k);
    hv_R_int = hv_L_int + V(3,k) * alpha(k);

    % Compute wave speeds and velocity of the intermediate state
    c_L_int = sqrt( g * h_L_int );
    c_R_int = sqrt( g * h_R_int );

    v_L_int = hv_L_int / h_L_int;
    v_R_int = hv_R_int / h_R_int;

    % Compute characteristic speeds (eigenvalues of the Jacobian)
    lambda_L_int(1) = v_L_int;
    lambda_L_int(2) = v_L_int + c_L_int;
    lambda_L_int(3) = v_L_int - c_L_int;

    lambda_R_int(1) = v_R_int;
    lambda_R_int(2) = v_R_int + c_R_int;
    lambda_R_int(3) = v_R_int - c_R_int;

    % Detect "transsonic rarefaction" (= problematic situation)
    if lambda_L_int(k) < 0 && lambda_R_int(k) > 0
        numerator = ( lambda_L_int(k) + lambda_R_int(k) ) * lambda(k) ...
            - 2 * lambda_L_int(k) * lambda_R_int(k);
        denominator = lambda_R_int(k) - lambda_L_int(k);

        lambda_fix(k) = numerator / denominator;
    else
        lambda_fix(k) = abs(lambda(k));
    end
end

end
