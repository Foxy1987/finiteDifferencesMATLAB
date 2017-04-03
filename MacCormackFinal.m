function [U]  = MacCormackFinal(L, nx)

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
    
    % set up vectors to store flow-field variables at throat of nozzle
    T = tempOld(16);
    RHO = rhoOld(16);
    M = vOld(16)/sqrt(tempOld(16));
    P = rhoOld(16)*tempOld(16);
    DRHOAVG = 0;
    DVAVG = 0;
    
    rhoPredicted =rhoOld;
    vPredicted = vOld;
    tempPredicted = tempOld;

    % since dt varies depending on the speed of the wave, we calculate
    % within the loop
    courantNum = 0.5;
    gamma = 1.4;
    i = 1;
    
    while (i < 3)
        % carry out MacCormack method for interior grid points 2:nx   
        
        dRho = -rhoOld(2:nx).*((vOld(3:nx1)-vOld(2:nx))./hx)-(rhoOld(2:nx).*vOld(2:nx)).*((log(A(3:nx1))-log(A(2:nx)))./hx)-vOld(2:nx).*((rhoOld(3:nx1)-rhoOld(2:nx))/hx);
        dV = -vOld(2:nx).*((vOld(3:nx1)-vOld(2:nx))./hx)-(1/gamma).*(((tempOld(3:nx1)-tempOld(2:nx))./hx)+          (tempOld(2:nx)/rhoOld(2:nx)).*((rhoOld(3:nx1)-rhoOld(2:nx))/hx));
        dTemp = -vOld(2:nx).*((tempOld(3:nx1)-tempOld(2:nx))/hx) - (gamma-1).*tempOld(2:nx).*(((vOld(3:nx1)-vOld(2:nx))/hx)+vOld(2:nx).*((log(A(3:nx1))-log(A(2:nx)))/hx));
        
        % find ht for all grid points and then choose the minimum as the
        % time step used in the MacCormack method
        % Calculate time step for ALL grid point on the current time step
        a = sqrt(tempOld);
        htVec = courantNum*(hx./(a + vOld));
        ht = min(htVec);
        
        % obtain predicted values 
        % we also have to deal with boundaries for predicted values
        
        rhoPredicted(2:nx) = rhoOld(2:nx) + dRho.*ht;
        vPredicted(2:nx) = vOld(2:nx) + dV.*ht;
        tempPredicted(2:nx) = tempOld(2:nx) + dTemp.*ht;
        
        % everything is correct up to here!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Predictor step is complete  %%%%%%%
        %%%%%%%%  Start the corrector step   %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Corrector step: backward difference equations 7.57-7.59
 
    
        dRhoBar = -rhoPredicted(2:nx).*((vPredicted(2:nx)-vPredicted(1:nx-1))/(hx)) - (rhoPredicted(2:nx)).*(vPredicted(2:nx)).*((log(A(2:nx))-log(A(1:nx-1)))./hx)-     vPredicted(2:nx).*((rhoPredicted(2:nx)-rhoPredicted(1:nx-1))./hx);
        dvBar = -vPredicted(2:nx).*((vPredicted(2:nx)-vPredicted(1:nx-1))./hx) - (1/gamma).*(((tempPredicted(2:nx)-tempPredicted(1:nx-1))./hx) + (tempPredicted(2:nx)./rhoPredicted(2:nx)).*((rhoPredicted(2:nx)-rhoPredicted(1:nx-1))./hx));
        dTempBar = -vPredicted(2:nx).*((tempPredicted(2:nx)-tempPredicted(1:nx-1))./hx) - (gamma-1).*tempPredicted(2:nx).*(((vPredicted(2:nx)-vPredicted(1:nx-1))/hx) + vPredicted(2:nx).*((log(A(2:nx))-log(A(1:nx-1)))/hx));
        
        % take average of the backward and forward difference using equations
        % 7.60-7.62
        dRhoAvg = (dRho + dRhoBar)./2;
        dvAvg = (dV + dvBar)./2;
        dTempAvg = (dTemp + dTempBar)./2;
        
        % obtain corrected values for p, V, T at time t + dt using equations
        % 7.63-7.65
       
        rhoNew(2:nx) = rhoOld(2:nx) + dRhoAvg.*ht;
        vNew(2:nx) = vOld(2:nx) + dvAvg.*ht;
        tempNew(2:nx) = tempOld(2:nx) + dTempAvg.*ht;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Corrector step is complete  %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
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
        T = [T tempNew(16)];
        M = [M vNew(16)/sqrt(tempNew(16))];   
        RHO = [RHO rhoNew(16)];
        P = [P pnew(16)];
        DRHOAVG = [DRHOAVG dRhoAvg(16)];
        DVAVG = [DVAVG dvAvg(16)];
        
        % set uOld = uNew;
        vOld = vNew;
        rhoOld = rhoNew;
        tempOld = tempNew;

        i = i + 1;
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
    
    U = [xvec', A', machNo', pnew', rhoNew', tempNew'];

end