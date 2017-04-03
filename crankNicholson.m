function [ U ] = crankNicholson( a, g1, g2, L, T, nx, nt )
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

    rvec = zeros(1, nx);

    % Boundary Condition at first and last column
    U(1, 1) = feval(g1, ht); %BC at t_(n+1)
    U(1, nx1) = feval(g2, ht);

    % call Stencil
    avec(2:nx) = -r/2;
    bvec(1:nx) = 1+r;
    cvec(1:nx-1) = -r/2;
 
    % FOR ALL OTHER TIME LEVELS
    for i = 2:nt1
        U(i, 1) = feval(g1, (i-1)*ht); 
        U(i, nx1) = feval(g2, (i-1)*ht);
        
        rvec(1) = (r/2)*U(i-1, 1) + (1-r)*U(i-1, 2) + (r/2)*U(i-1, 3) +(r/2)*U(i, 1);
        rvec(2:nx-1) = (r/2)*U(i-1, 1:nx-2) + (1-r)*U(i-1, 2:nx-1) + (r/2)*U(i-1, 3:nx);
        rvec(nx) = (r/2)*U(i-1, nx-1) + (1-r)*U(i-1, nx) + (r/2)*U(i, nx1);
        
        rvec = tridiag( avec, bvec, cvec, rvec);
        U(i, 2:nx) = rvec(2:nx);
    end

    surf(xvec, tvec, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
    axis tight

    camlight left
end

