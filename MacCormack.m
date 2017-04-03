function U = MacCormack(a, L, T, nx, nt)
  
nx1 = nx+1;
hx = L/nx;
nt1 = nt+1;
ht = T/nt;

R = (a*ht)/(hx);

xvec = hx*(0:nx);

tvec = ht*(0:nt);

% initial condition
U = zeros(nt1, nx1);
U(1,:) = (sin(pi*xvec(:)).^40);

uTemp = zeros(1, nx1);

for i = 1:nt
    % predictor step
    uTemp(2:nx) = U(i,2:nx) - R*(U(i, 3:nx1)-U(i,2:nx));
    % corrector step
    U(i+1, 2:nx) = (1/2)*((U(i, 2:nx) + uTemp(2:nx)) + R*(uTemp(2:nx) - uTemp(1:nx-1)));
end

figure
 [X, Y] = meshgrid(ht*(0:nt), hx*(0:nx));
 surf(X, Y, U','FaceColor','interp',...
    'EdgeColor','blue',...
    'FaceLighting','phong')
 axis tight
 camlight left
% 
 xlabel('space');
 ylabel('time');