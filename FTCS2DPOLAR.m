function [ U ] = FTCS2DPOLAR( R, TH, T, Mr, Mt, nt )

    %set grid dimensions
    Mt=floor(Mt);
    Mr1 = Mr+1;
    hr = R/Mr;
    Mt1 = floor(Mt+1);
    htheta = TH/Mt;
    nt1 = nt+1;
    ht = T/nt;
    
    % create r, theta, t vectors 
    rVec = hr*(0:Mr);
    thetaVec = htheta*(0:Mt);
    tVec=ht*(0:nt);
    
    % create solution vector
    U = zeros(Mr1, floor(Mt1), nt1); % cnx+1 because we run from 0 to nx not 1 to nx
    
    % set initial conditions
    U(:,:, 1) =0;
    
    % set boundary conditions all set to 0
   
    U(1,:, :) = sin(4*thetaVec')*sin(tVec);
    
    
    for jt = 1:nt-1
        % loop r find 2nd derivative of theta
        for j = 2:Mr
            % fix theta find first derivative of r
            for k = 2:Mt
                U(j, k, jt+1) = U(j, k, jt) + (1/rVec(j))*(ht/hr^2)*(((rVec(j)+rVec(j+1))/2)*(U(j+1, k, jt) ... 
                -U(j, k, jt)) + ((rVec(j)+rVec(j-1))/2)*(U(j, k, jt) -U(j-1, k, jt))) + ... 
                (1/rVec(j)^2)*(ht/htheta^2)*(U(j, k+1, jt) -2*U(j, k, jt) + U(j, k-1, jt));
            end
            U(j, 1, jt+1) = U(j, 1, jt) + (1/rVec(j))*(ht/hr^2)*(((rVec(j)+rVec(j+1))/2)*(U(j+1, 1, jt) ... 
                -U(j, 1, jt)) - ((rVec(j)+rVec(j-1))/2)*(U(j, 1, jt) -U(j-1, 1, jt))) + ... 
                (1/rVec(j)^2)*(ht/htheta^2)*(U(j, 2, jt) -2*U(j, 1, jt) + U(j, Mt, jt));
        end

        U(1, 1:Mt, jt+1) = (1-(4*ht)/hr^2)*U(1, 1:Mt, jt) + (2*htheta*ht/pi*hr)*sum(U(1, 1:Mt,jt));
    end

    
     mov(1:100/2) = struct('cdata', [],...
                'colormap', []);
     figure
     axis([0 1 0 1 -10 10]);
     [X,Y] = meshgrid(0:hr:1, 0:htheta:2*pi);

      for jt = 1:100  
        surf(X, Y, U(:, :,jt)','FaceColor','interp',...
       'EdgeColor','blue',...
       'FaceLighting','phong'); camlight('left');
         colormap (jet(64));
         colorbar('vertical');
         axis([0 1 0 1 -100 100]);
         drawnow; pause(0.001);
 
         mov(jt) = getframe(gcf);
     end
      
    % movie2avi(mov,'so1.avi','compression', 'None','fps', 10000);
 
    % provide 2 dimensional slices of the solution (holding either x or y constant while the other varies)