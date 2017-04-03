function [ U ] = implicit2( a, L, T, nx, nt )
    
    nx1 = nx+1;
    hx = L/nx;
    nt1 = nt+1;
    ht = T/nt;

    R = (a*ht)/(2*hx);
    xvec = hx*(0:nx);
    tvec=ht*(0:nt);

    % initial condition
    U = zeros(nt1, nx1); % cnx+1 because we run from 0 to nx not 1 to nx
    U(1,:) = sin(pi*xvec(:)).^80;
    %feval(u0, xvec(:)');
    
    % form circulant matrix with periodic boundary conditions
    % form Q matrix (This is tridiagonal so just store a, b, c vectors)
    B_avec(2:nx-1) = -R/2;
    B_bvec(1:nx-1) = 1;
    B_cvec(1:nx-2) = R/2;
    %B= diag(ones(nx-1, 1), 0) + diag(R/2.*ones(nx-2, 1), 1) +diag(-R/2.*ones(nx-2, 1), -1);
    % change appropriate elements to get tridiagonal matrix B
    B_bvec([1 nx-1]) = [1+(R/2) 1-(R/2)];
    
    % make circulant matrix tridiagonal with Sherman-Morrison algorithm   
    zVec = zeros(nx-1, 1); wVec = zeros(nx-1, 1);
    zVec([1 nx-1], 1) = R/2; wVec([1 nx-1], 1) = [1 -1];
   
    for i = 1:nt  
        
        U(i, 1) = U(i, nx);
        % RHS is just U_k
        rvec = [U(i,2:nx)' wVec];
        
        % We pass 2 into tridiag2 because there are 2 RHS vectors
        rvec =tridiag2(B_avec, B_bvec, B_cvec, rvec, 2);
        %y1 = B\rvec(:, 1);
        %y2 = B\rvec(:, 2);
        
        
        % now y1 and y2 are stored in r. Use them to find beta
        y1 = rvec(:,1); y2 = rvec(:, 2);
        beta = (zVec.*y1)./(1-(zVec.*y2));
    
        U(i+1, 2:nx) = y1 + y2.*beta;
        
        
    end
   
    surf(xvec, tvec, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
    axis tight

    camlight left
end
