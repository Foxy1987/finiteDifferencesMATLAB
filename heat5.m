% question 6
function U = heat5(a, f, u0, g1, g2, L, T, nx, nt)
  
% parameters
% ghost point adds 1 more point so we have N+2 equations
nx2 = nx+2;
hx = L/nx;
nt1 = nt+1;
ht = T/nt;

r = (a*ht)/(hx^2);
r1 = 1- 2*r; 

xvec = hx*(0:nx);
%xvec = 0:hx:L+hx;

tvec = ht*(0:nt);
%tvec = 0:ht:T+ht

% initial condition
U = zeros(nt1, nx2);
%U(1,1:nx+1) = feval(u0, xvec(:)');

for k = 1:nx+1
   U(1, k) = sin(2*pi*xvec(k)); 
end

U(:, 1) = sin(2*pi*tvec);

% Boundary Condition at first and last column
%U(2:nt1, 1) = feval(g1, (1:nt)*ht);
%U(2:nt1, nx2) = feval(g2, (1:nt)*ht);

U(1:nt1, nx2) = U(1, nx) + 4*pi*hx;

for k = 1:nt
    % from the column adjacent to the rightmost boundary all the way to the
    % second column
    for j = 2:nx
        U(k+1,j) = U(k, j) + r*(U(k, j+1) -2*U(k, j) + U(k, j-1));
    end
    
    U(k+1, j+2) = U(k+1, j) + 4*pi*hx;
end


%  xx = linspace(0, 1, 22)
%  plot(xx, U(3, :))
% 
% 
 figure
 [X, Y] = meshgrid(ht*(0:nt), hx*(0:nx+1));
 surf(X, Y, U','FaceColor','interp',...
    'EdgeColor','blue',...
    'FaceLighting','phong')
 axis tight
 camlight left
% 
 xlabel('space');
 ylabel('time');
 zlabel('Temperature');


end