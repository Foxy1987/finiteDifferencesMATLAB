function [ U ] = btcs( a, g1, g2, L, T, nx, nt )
% define size of vectors nx-1 + BC = nx+1
    nx1 = nx+1;
    hx = L/nx;
    nt1 = nt+1;
    ht = T/nt;

    r = (a*ht)/(hx^2);

    xvec = hx*(0:nx);
    tvec=ht*(0:nt);

    % initial condition
    U = zeros(nt1, nx1);
    U(1,:) = sin(2*pi*xvec(:));
    %feval(u0, xvec(:)');

    rvec = zeros(1, nx-1);

    % Boundary Condition at first and last column
    U(2, 1) = feval(g1, ht); %BC at t_(n+1)
    U(2, nx1) = feval(g2, ht);
    
    % compute the solution at t = ht: use the program tridiag to solve the
    % system Qu = r at t = ht 
    % AT TIME LEVEL 1

    % call RHS with Boundary conditions (if BC set to 0 ignore)
    rvec(1) = U(1,2) + r*U(2,1);
    rvec(2:nx-2) = U(1, 3:nx-1);
    rvec(nx-1) = U(1, nx) + r*U(2, nx1);


    % call Stencil
    avec(2:nx-1) = -r;
    bvec(1:nx-1) = 1+2*r;
    cvec(1:nx-2) = -r;

    % call tridiag
    rvec = tridiag( avec, bvec, cvec, rvec);
    U(2, 2:nx) = rvec(1:nx-1);
    

    % FOR ALL OTHER TIME LEVELS
    for k = 2:nt1
            % calculate U at the boundaries
            U(k, 1) = feval(g1, (k-1)*ht); %BC at t_(n+1)
            U(k, nx1) = feval(g2, (k-1)*ht);  
            
            % calculate the RHS
            rvec(1) = U(k-1,2) + r*U(k,1);
            rvec(2:nx-2) = U(k-1, 3:nx-1);
            rvec(nx-1) = U(k-1, nx) + r*U(k, nx1);
            
            % find the solution vector u at time t_{n+1} by solving Qu = r
            rvec = tridiag( avec, bvec, cvec, rvec);
            U(k, 2:nx) = rvec(1:nx-1);
    end
    
    surf(xvec, tvec, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
    axis tight

    camlight left
  
end