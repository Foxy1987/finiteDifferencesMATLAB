% question 6
function [U, tvec, xvec] = LeapFrogWave(a, L, T, nx, nt)
  
    % parameters
    % ghost point adds 1 more point so we have N+2 equations
    nx1 = nx+1;
    hx = L/nx;
    ht = T/nt;
    nt1 = nt+1;
    R = (a*ht)/(hx);
    xvec = hx*(0:nx);
    tvec = ht*(0:nt);
    
    % initial condition
    U = zeros(nt1, nx1);
    for k = 2:nx
        U(1, k) = sin(pi*xvec(k)).^80; 
    end
    
     % boundary conditions
    U(:, 1) = U(:, nx1);

    % kick the leap frog off with a one-step explicit method
    % Use Lax-Wendroff method since it is O(x^2) + O(t^2)
    % for lax-wendroff make r^2 >= R^2/4
    U(2, 2:nx) = U(1, 2:nx) - (R/2).*(U(1, 3:nx1) - U(1:nx-1)) + ((R^2)/2).*(U(1, 1:nx-1) ...
                            - 2.*U(1, 2:nx) + U(1, 3:nx1));
    
    for k = 2:nt
        U(:, 1) = U(:, nx1);
    	U(k+1, 2:nx) = U(k-1, 2:nx) - R.*(U(k, 3:nx1) - U(k, 1:nx-1));  
    end


     figure;
     [X, Y] = meshgrid(ht*(0:nt), hx*(0:nx));
     surf(Y, X, U','FaceColor','interp',...
        'EdgeColor','blue',...
        'FaceLighting','phong')
     axis tight
     camlight left
     
     xlabel('space');
     ylabel('time');
     zlabel('Temperature');


end