function [tout1, tout2, Nout1, Nout2] = RK4_test(tsim, DT)
%RK4_TEST Summary of this function goes here
%   Detailed explanation goes here

    N0 = zeros(3,1);

    % Time span (yn = 1/tau, where tau is the lifetime)
    maxt = tsim;
    mint = 0.0;
    npts = 4001; 
    dt = (maxt - mint)/(npts - 1.0);

    tspan = mint:dt:maxt;
    
    param = 1;

    odefun = @(t, N) test(t, N, param); % Anonymous handle to function

    reltol = 1E-6;
    options = odeset('RelTol', reltol);

    % Runge-Kutta implementation
    [tout1, Nout1] = ode45(odefun, tspan, N0, options); 
    
    figure
    hold on
    title('ode45')
    plot(tout1, Nout1(:,1))
    plot(tout1, Nout1(:,2))
    plot(tout1, Nout1(:,3))
    grid on
    
    [tout2, Nout2] = RK4(odefun, N0, mint, maxt, DT);
    
    figure
    hold on
    title('RK4')
    plot(tout2, Nout2(:,1))
    plot(tout2, Nout2(:,2))
    plot(tout2, Nout2(:,3))
    grid on
    
    return

    function dN = test(t, N0, ~)
           
        %y1 = N0(1);
        %y2 = N0(2);
        %y3 = N0(3);
        
        dN = zeros(size(N0));
        
        % dy1/dt
        dN(1) = sin(t);
        
        % dy2/dt
        dN(2) = cos(t);
        
        % dy3/dt
        dN(3) = 0;
                
    end
end

