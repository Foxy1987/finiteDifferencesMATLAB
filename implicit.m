%implicit method for solving the viscous burgers' equation
function [ U ] = implicit( a, L, T, nx, nt )
    
    nx1 = nx+1;
    hx = L/nx;
    nt1 = nt+1;
    ht = T/nt;

    r = (a*ht)/(hx^2);
    R = ht/(2*hx);
    
    xvec = hx*(0:nx);
    tvec=ht*(0:nt);

    % initial condition
    U = zeros(nt1, nx1); % cnx+1 because we run from 0 to nx not 1 to nx
    U(1,:) = sin(2*pi*xvec(:));
    %feval(u0, xvec(:)');

    fvec = zeros(1, nx-1);
    
    % Boundary Condition at first and last column
    U(:, 1) = sin(2*pi*tvec(:));%zeros(nt1, 1); % boundary at k = 0 or in matlab k = 1
    U(:, nx1) = zeros(nt1, 1); % boundary at k = M or in matlab k=M+1 (or nx1)
 
    % march through time levels
    for i = 1:nt  
        
        % set uOld at beginning of time step
        error = inf; relTol = 1e-5; maxIterates = 30;
        itCount = 0;
        % initial condition and then previous time step is our initial guess for the Newton method
        uOld=U(i, :);
        uNew0=U(i, :);
        % solve nx-1 nonlinear equations in nx-1 unknowns f(x) = 0
        % The diagonal entries of the Jacobian are given in a, b, c
        % Only output uNew when the error is less than 1e-5 or we reach max
        % number of iterates
        while error >= relTol && itCount < maxIterates
            
            avec(2:nx-1)=(-R/2)*uOld(2:nx-1) - r;
            bvec(1:nx-1)= 1+(R/2)*(uOld(3:nx+1) - uOld(1:nx-1)) + 2*r;
            cvec(1:nx-2) = (R/2)*uOld(1:nx-2) - r;
          
            % FIRST THREE TERMS COME FROM THE LAST NEWTON ITERATE
            % LAST TERM WILL ALWAYS BE UOLD
            fvec(1:nx-1) = uNew0(2:nx) + (R/2)*(uNew0(3:nx+1) - uNew0(1:nx-1)) - r*(uNew0(1:nx-1) - 2*uNew0(2:nx) + uNew0(3:nx+1)) - uOld(2:nx);
            
            % fvec(uOld) and [avec; bvec; cvec] = f and f' respectively
            % F'(uOld)*deltaU = -F(uOld) => solve for deltaU
            % We solve F'(uOld)deltaU = -F(uOld) (which is Ax = b) using the Thomas Algorithm
            % since F' is tridiagonal
            deltaU = tridiag(avec, bvec, cvec, fvec);
            
            uNew1(2:nx) =  uNew0(2:nx) - deltaU;
            error = norm(deltaU, inf);
            uNew0(2:nx) = uNew1(2:nx);
            
            if itCount == maxIterates
                disp(' ')
                disp('*** The answer may possibly not satisfy the correct tolerance ****');
            end
            
        end
        U(i+1, 2:nx) = uNew1(2:nx);
       
    end

    surf(xvec, tvec, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
    axis tight

    camlight left
end
