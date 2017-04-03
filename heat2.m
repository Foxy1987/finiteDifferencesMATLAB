% second order leapfrog question 2
function U = heat2(a, f, u0, g1, g2, L, T, nx, nt)
  
% parameters
% ghost point adds 1 more point so we have N+2 equations
nx2 = nx+2;
hx = L/nx;
nt1 = nt+1;
ht = T/nt;

r = (a*ht)/(hx^2);
r1 = 1- 2*r; 

xvec = hx*(0:nx+1);

tvec = ht*(0:nt);

% initial condition
U = zeros(nt1, nx2);
%U(1,:) = feval(u0, xvec(:)');

% the one cell on the end is the ghost cell
for i = 1:nx+2
   U(1, i) = cos(pi*xvec(i)/2); 
end

% Boundary Condition at first and last column
%U(2:nt1, 1) = feval(g1, (1:nt)*ht);
U(1:nt1, nx2) = 0;%feval(g2, (1:nt)*ht);

% In leapfrog current solution requires solutions at 2 previous times
% we already have 1 solution for all x and 1st time
% get solution for all x for second time
for j = 1:nx2 
    U(2, j) = cos((pi*(j-1)*hx)/2)*exp(-a*((pi^2)*(2)*ht)/4);
end

for k = 2:nt
    % from the column adjacent to the rightmost boundary all the way to the
    % second column
    for j = nx+1:-1:2
        U(k+1,j) = U(k-1, j) + r*(U(k, j+1) -2*U(k, j) + U(k, j-1));
    end
    % only cell left to fill in the row is 1 and 3
    % when we get to the second column set U(k, -1) = U(k, 1)
    U(k+1, 1) = U(k+1, 3);
end


figure
[X, Y] = meshgrid(hx*(0:nx+1), ht*(0:nt));
surf(X, Y, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
axis tight
%camlight left

xlabel('space');
ylabel('time');
zlabel('Temperature');
