function U = MacCormackFinal2(L, nx)

    format long
    % N+1 points 2:N points are the interior grid points and grid points 1 and
    % N+1 are the boundary layers
    nx1 = nx+1;
    hx = L/nx;
    xvec = hx*(0:nx); % This gives us N+1 points because because MATLAB starts at position 1 not 0

    % compute nozzle shape
    A = 1+2.2.*(xvec - 1.5).^2;
    
    % set initial conditions for p, V, T
    rhoOld = 1-0.3146*xvec;
    tempOld = 1-0.2314*xvec;
    vOld = (0.1+1.09.*xvec).*sqrt(tempOld);
    rhoNew = rhoOld; vNew = vOld; tempNew=tempOld;
    
    dRho=zeros(1, length(2:nx)); dV=zeros(1, length(2:nx)); dTemp=zeros(1, length(2:nx));
    
    dRhoBar=zeros(1, length(2:nx)); dvBar=zeros(1, length(2:nx)); dTempBar=zeros(1, length(2:nx));
    rhoPredicted =rhoOld;
    vPredicted = vOld;
    tempPredicted = tempOld;
    
    % set up vectors to store flow-field variables at throat of nozzle
    T = zeros(1, 1400);
    T(1) = tempOld(15);
    RHO = rhoOld(15);
    M = vOld(15)/sqrt(tempOld(15));
    P = rhoOld(15)*tempOld(15);
    DRHOAVG = 0;
    DVAVG = 0;

    % since dt varies depending on the speed of the wave, we calculate
    % within the loop
    courantNum = 0.5;
    gamma = 1.4;
    i = 1;
    while (i < 1400)
        i = i + 1;
        % carry out MacCormack method for interior grid points 2:nx   
        for j =2:nx
            dRho(j-1) = -rhoOld(j).*((vOld(j+1)-vOld(j))./hx)-(rhoOld(j).*vOld(j)).*((log(A(j+1))-log(A(j)))./hx)-vOld(j).*((rhoOld(j+1)-rhoOld(j))/hx);
            dV(j-1) = -vOld(j).*((vOld(j+1)-vOld(j))./hx)-(1/gamma).*(((tempOld(j+1)-tempOld(j))./hx)+(tempOld(j)/rhoOld(j)).*((rhoOld(j+1)-rhoOld(j))/hx));
            dTemp(j-1) = -vOld(j).*((tempOld(j+1)-tempOld(j))/hx) - (gamma-1).*tempOld(j).*(((vOld(j+1)-vOld(j))/hx)+vOld(j).*((log(A(j+1))-log(A(j)))/hx));
        end
        
        a = sqrt(tempOld);
        htVec = courantNum*(hx./(a + vOld));
        ht = min(htVec);
        
        rhoPredicted(2:nx) = rhoOld(2:nx) + dRho.*ht;
        vPredicted(2:nx) = vOld(2:nx) + dV.*ht;
        tempPredicted(2:nx) = tempOld(2:nx) + dTemp.*ht;

      
        for j = 2:nx
            dRhoBar(j-1) = -rhoPredicted(j).*((vPredicted(j)-vPredicted(j-1))/(hx)) - (rhoPredicted(j)).*(vPredicted(j)).*((log(A(j))-log(A(j-1)))./hx)-     vPredicted(j).*((rhoPredicted(j)-rhoPredicted(j-1))./hx);
            dvBar(j-1) = -vPredicted(j).*((vPredicted(j)-vPredicted(j-1))./hx) - (1/gamma).*(((tempPredicted(j)-tempPredicted(j-1))./hx) + (tempPredicted(j)./rhoPredicted(j)).*((rhoPredicted(j)-rhoPredicted(j-1))./hx));
            dTempBar(j-1) = -vPredicted(j).*((tempPredicted(j)-tempPredicted(j-1))./hx) - (gamma-1).*tempPredicted(j).*(((vPredicted(j)-vPredicted(j-1))/hx) + vPredicted(j).*((log(A(j))-log(A(j-1)))/hx));
        end

       
        dRhoAvg = (dRho + dRhoBar)./2;
        dvAvg = (dV + dvBar)./2;
        dTempAvg = (dTemp + dTempBar)./2;


        rhoNew(2:nx) = rhoOld(2:nx) + dRhoAvg.*ht;
        vNew(2:nx) = vOld(2:nx) + dvAvg.*ht;
        tempNew(2:nx) = tempOld(2:nx) + dTempAvg.*ht;

            
        % Now we have to calculate the flow-field variables at the boundary
        % points. 
        % extrapolate at the subsonic inflow boundary
        vNew(1) = 2*vNew(2)-vNew(3);
        rhoNew(1)=1; tempNew(1)=1;
        % extrapolate at the supersonic outflow boundary 
        rhoNew(nx1) = 2*rhoNew(nx) - rhoNew(nx-1);
        vNew(nx1) = 2*vNew(nx) - vNew(nx-1);
        tempNew(nx1) = 2*tempNew(nx) - tempNew(nx-1);
        
        % define a nondimensional pressure as the local static pressure
        % divided by the reservoir pressure p0, the equation of the state
        % is given by: p=rho*T. p at grid point 16 should = 0.349
        pnew = rhoNew.*tempNew;
        machNo = vNew./sqrt(tempNew);
        
        % append throat values to vectors so we can shock
        T(i) = tempNew(15);
        M = [M vNew(15)/sqrt(tempNew(15))];   
        RHO = [RHO rhoNew(15)];
        P = [P pnew(15)];
        DRHOAVG = [DRHOAVG dRhoAvg(15)];
        DVAVG = [DVAVG dvAvg(15)];
        
        % set uOld = uNew;
        vOld = vNew;
        rhoOld = rhoNew;
        tempOld = tempNew;

        
    end
    
    % plot timewise variations of flow-field variables at throat 
    figure
    subplot(2, 2, 1);
    plot(T, '-b', 'LineWidth', 4)
    ylabel('T/T0');
    subplot(2, 2, 2);
    plot(M, '-r', 'LineWidth', 4)
    ylabel('M')
    subplot(2, 2, 3);
    plot(RHO, '-c', 'LineWidth', 4)
    ylabel('rho/rho0')
    subplot(2, 2, 4);
    plot(P, '-g', 'LineWidth', 4)
    ylabel('p/p0');
    
    % plot timewise variations in the time-derivates of density and velocity
    figure;
    semilogy(abs(DRHOAVG'), '-g', 'LineWidth', 4)
    hold on
    semilogy(abs(DVAVG'), '-r', 'LineWidth', 4)
    legend('dRhoAvg', 'dvAvg');
    
    %  plot steady-state flow-field variables along the length of the nozzle 
    figure
    subplot(2, 2, 1);
    plot(xvec, machNo);
    subplot(2, 2, 2);
    plot(xvec, pnew) 
    subplot(2, 2, 3);
    plot(xvec, rhoNew)
    subplot(2, 2, 4)
    plot(xvec, tempNew)

end