function [t, Nout] = RK4(func, Nin, t0, t1, dt, varargin)
%RK4 Implementation of a 4th order Runge Kutta solver for a fixed time step
%
%Syntax:
%
%	[Nout, t] = RK4(func, Nin, t0, t1, dt, [varargin]);
%
%Arguments:
%
%	func        Function handle (passed via func = @somefunction)
%   Nin         Array of initial values (npts X nvar)
%   t0          Initial time
%   t1          Final time
%   dt          Time step
%
%Return values:
%
%   Nout        Solution array
%   t           Time array

    % Test if func is a function handle
    if (not(isa(func, 'function_handle')))
   
        error('First argument to RK4 is not a function handle')

    end
    
    % Get dimensions of input array
    s = size(Nin);
    
    if (length(s) > 2)
        
        error('Input array must be one dimensional')
        
    end
    
    if (and(not(s(1)==1),not(s(2)==1)))
        
        error('Input array must be one dimensional')
        
    end
    
    if (s(1)==1)
        
        nvar = s(2);
        
    else
        
        nvar = s(1);
        % Transpose array
        Nin = transpose(Nin);
        
        
    end
    
    % Set simulation time
    tsim = t1 - t0;
    
    % Number of points to calculate
    npts = ceil(abs(tsim/dt)) + 1;
    
    % Recalculate dt to make simulation time an integral multiple
    dt = tsim/(npts - 1);
    
    % Set up time array
    t = t0:dt:t1;
    
    % Transpose time array
    t = transpose(t);
    
    % Set up output array
    Nout = zeros(npts, nvar);
    
    % Set initial values
    Nout(1,:) = Nin;
    
    for n = 1:npts-1
        
        % Calculate next step via RK4
        N0 = Nout(n,:);
        t0 = t(n);
        k1 = func(t0, N0)*dt; 
        k2 = func(t0 + 0.5*dt, N0 + 0.5*k1)*dt;
        k3 = func(t0 + 0.5*dt, N0 + 0.5*k2)*dt;
        k4 = func(t0 + dt, N0 + k3)*dt;
        Nout(n+1,:) = N0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
        
    end
    

end

