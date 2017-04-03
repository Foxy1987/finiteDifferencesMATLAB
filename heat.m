function U = heat(a, f, g1, g2, L, T, nx, nt)

% forward difference scheme for solving the initial boundary value problem
% of the heat equation


% U: the solution (u_{k, i}), k the time grid point index, i the space grid
% point index


% parameters
nx1 = nx+1;
hx = L/nx;
nt1 = nt+1;
ht = T/nt;

r = a*ht/(hx^2);
r1 = 1- 2*r; 

xvec = hx*(0:nx);
tvec = ht*(0:nt);

% initial condition
U = zeros(nt1, nx1);
U(1,:) = sin(2*pi*xvec(:));

% advance time
% boundary conditions
U(2:nt1, 1) = feval(g1, (1:nt)*ht);
U(2:nt1, nx1) = feval(g2, (1:nt)*ht);

for k = 1:nt
   U(k+1, 2:nx) = r*(U(k, 1:nx-1)+U(k, 3:nx+1))+r1*U(k, 2:nx) ...
       +ht*feval(f, xvec(2:nx), tvec(k));
end

figure
surf(xvec, tvec, U,'FaceColor','interp',...
   'EdgeColor','blue',...
   'FaceLighting','phong')
axis tight
view(-50,50)
camlight left

xlabel('space');
ylabel('time');
zlabel('Temperature');