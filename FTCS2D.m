function [ xvec, yvec, U ] = FTCS2D( a, L, W, T, nx, ny, nt )
    
    nx1 = nx+1;
    hx = L/nx;
    ny1 = ny+1;
    hy = W/ny;
    nt1 = nt+1;
    ht = T/nt;
    
    rx = a*ht/hx^2;
    ry = a*ht/hy^2;
    
    xvec = hx*(0:nx);
    yvec = hy*(0:ny);
    tvec=ht*(0:nt);
    
    U = zeros(nx1, ny1, nt1); % cnx+1 because we run from 0 to nx not 1 to nx
    
    % set initial conditions
    U(:,:, 1) = sin(pi*xvec')*cos(2*pi*yvec);  
    % set boundary conditions all set to 0
    U(:, 1, :) = 0;
    U(:, ny1, :) = 0; 
    U(1,:, :) = 0;
    U(nx1,:,:) = 0;
    dxx = zeros(nx1, ny1);
    dyy = zeros(nx1, ny1);
    
    for jt = 1:nt  
        % starting from y boundary
        % keep y fixed and calculate dxx
        for jy = 2:ny
            dxx(2:nx, jy) = rx*(U(1:nx-1, jy, jt)-2*U(2:nx, jy, jt)+U(3:nx+1, jy, jt));
        end
        %keep x fixed and calculate dyy
        for jx = 2:nx
            dyy(jx, 2:ny) = ry*(U(jx, 1:ny-1, jt)-2*U(jx, 2:ny, jt)+U(jx, 3:ny+1, jt));
        end
        
        U(2:nx, 2:ny, jt+1) = U(2:nx, 2:ny, jt) + (dyy(2:nx, 2:ny) + dxx(2:nx, 2:ny));
        
    end
    mov(1:10000/2) = struct('cdata', [],...
                        'colormap', []);
    figure
    axis([0 1 0 1 -1 1]);
    [X,Y] = meshgrid(0:.05:3, 0:.05:3);
    for jt = 1:10000  
        surf(X, Y, U(:, :, jt),'FaceColor','interp',...
       'EdgeColor','blue',...
       'FaceLighting','phong'); camlight('left');
        colormap (jet(64));
        colorbar('vertical');
        axis([0 3 0 3 -1 1]);
        drawnow; pause(0.01);

        mov(jt) = getframe(gcf);
    end
     
     movie2avi(mov,'so1.avi','compression', 'None','fps', 10000);
     
     
    % provide 2 dimensional slices of the solution (holding either x or y constant while the other varies)