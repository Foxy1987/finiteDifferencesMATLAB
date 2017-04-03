% question 6
function U = burgers(a, L, T, nx, nt)
  
    % parameters
    % ghost point adds 1 more point so we have N+2 equations
    nx1 = nx+1;
    hx = L/nx;
    ht = T/nt;

    r = (a*ht)/(hx^2);
    s = ht/(2*hx);
    r1 = 1- 2*r; 

    xvec = hx*(0:nx);

    % initial condition
    U = zeros(nt, nx1);
    for k = 1:nx+1
        U(1, k) = sin(2*pi*xvec(k)); 
    end
    
     % boundary conditions
    U(:, 1) = 0;%feval(g1, (1:nt)*ht);
    U(:, nx1) = 0;%feval(g2, (1:nt)*ht);

    for k = 1:nt
        for j = 2:nx
            U(k+1,j) = U(k, j) - s*((U(k, j+1) - U(k, j-1))*U(k, j)) + r*(U(k, j+1) - 2*U(k, j) + U(k, j-1));
        end
    end


     figure
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