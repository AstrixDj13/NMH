function [ V, lambda ] = compute_eigensystem_of_jacobian_x(c, u, v)
%COMPUTE_EIGENSYSTEM_OF_JACOBIAN_X Evaluate the eigenvectors and eigenvalues of
%the Jacobian matrix for the x-direction of the shallow water equations

% Allocate results
%   This is not ideal from a programming point of view since these arrays are
%   allocated and deallocated in every pass through the loop.
V = zeros(3,3);
lambda = zeros(3,1);

% Right eigenvectors
V(1,1) = 0.0;
V(2,1) = 0.0;
V(3,1) = 1.0;
V(1,2) = 1.0;
V(2,2) = u + c;
V(3,2) = v;
V(1,3) = 1.0;
V(2,3) = u - c;
V(3,3) = v;

% Eigenvalues (characteristic speeds)
lambda(1) = u;
lambda(2) = u + c;
lambda(3) = u - c;

end
