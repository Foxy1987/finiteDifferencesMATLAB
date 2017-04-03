function [ root ] = newtonSystem(u_init, maxIterates, relTol, fsys, a, b, c)
%% 
% u_init is the vector of initial guess
% maxIterates is the maximum number of iterations
% relTol is the error tolerance
%%

u0 = zeros(n, 1);
u0(:) = u_init(:);
error = inf;
itCount = 0;

while error >= relTol & itCount < maxIterates
    
    itCount = itCount + 1;
    rhs = fsys(u0);
    A = deriv_fsys(u0);
    % we can use the mldivide \ operator to solve the system of linear
    % equations Ax = B. A and B must have the same number of rows
    delta = A\rhs;
    
    % OR we can use Thomas algorithm to solve the tridiagonal system using
    % the just sub, main, and superdiagonal (a, b, c)
    % delta_U will be returned by the trid function
    
    u1 = u0 - delta;
    error = norm(delta, inf);
    
    u0 = u1;
end

% return the solution vector
root = u1;
