function [ V, lambda ] = compute_eigensystem_of_jacobian_y(c, u, v)
%COMPUTE_EIGENSYSTEM_OF_JACOBIAN_Y Evaluate the eigenvectors and eigenvalues of
%the Jacobian matrix for the y-direction of the shallow water equations

% Allocate results
%   This is not ideal from a programming point of view since these arrays are
%   allocated and deallocated in every pass through the loop.
V = zeros(3,3);
lambda = zeros(3,1);

% Right eigenvectors
V(1,1) = 0.0;
V(2,1) = 1.0;
V(3,1) = 0.0;
V(1,2) = 1.0;
V(2,2) = u;
V(3,2) = v + c;
V(1,3) = 1.0;
V(2,3) = u;
V(3,3) = v - c;

% Eigenvalues (characteristic speeds)
lambda(1) = v;
lambda(2) = v + c;
lambda(3) = v - c;

end
