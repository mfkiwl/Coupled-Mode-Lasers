function [tout, Nout] = runSimulation(tsim, N0, param, varargin)
%RUNSIMULATION Solves the rate equations using the Runge-Kutta method 
%   
%% Description
%
%% Usage
%
%   [tout, Nout] = runSimulation(tsim, N0, param, varargin);
%
%% Arguments
%
%	tsim        integer specifying the simulation time tspan via tspan =
%               tsim/yn, where yn is the recombination rate (in 1/ns) 
%
%	N0          Vector of inital conditions for the equations
%
%               Carrier concentrations:
%                   MA   = N0(1)        Carrier concentration in guide A
%                   MB   = N0(2)        Carrier concentration in guide B
%               Optical fields:
%                   YA   = N0(3)        Amplitude in guide A
%                   YB   = N0(4)        Amplitude in guide B
%                   phiA = N0(5)        Phase of optical field in guide A
%                   phiB = N0(6)        Phase of optical field in guide B
%
%
%   param       structure containing laser parameters
%
%               param.yn        carrier recombination rate (1/tau_N)
%               param.kp        cavity loss rate (1/(2*tau_p))
%               param.aH        Linewidth enhancement factor
%               param.kinj      Injection level rate (K_inj/tau_N)
%               param.DWA       Free-running minus laser A angular frequency 
%               param.DWB       Free-running minus laser B angular frequency 
%               param.DW        Cavity detuning (see notes)
%               param.DWinj 	Injection detuning (W_inj - W)
%               param.QA        Normalised pumping rate in guide A    
%               param.QB        Normalised pumping rate in guide A 
%               param.eta       Amplitude of coupling coefficient (not used)
%               param.etaAB     Amplitude of coupling coefficient AB (asymmetric)
%               param.etaBA     Amplitude of coupling coefficient BA (asymmetric)
%               param.theta     Phase of coupling coefficient 
%
%   varargin    optional values:
%
%
%% Returns
% 
%   tout        4001 x 1 array of time values
%
%   Nout        4001 x 5 array containing the time evolution of variables
%
%                   Nout(:,1)   Carrier concentration in guide A
%                   Nout(:,2)   Carrier concentration in guide B
%                   Nout(:,3)   Optical amplitude in guide A
%                   Nout(:,4)   Optical amplitude in guide B
%                   Nout(:,5)   Phase of optical field in guide A
%                   Nout(:,6)   Phase of optical field in guide B
%
%% Code

    % Time span (yn = 1/tau, where tau is the lifetime)
    yn = param.yn;
    maxt = tsim/yn;
    mint = 0.0;
    npts = 4001; 
    dt = (maxt - mint)/(npts - 1.0);

    tspan = mint:dt:maxt;

    odefun = @(t, N) coupledInj(t, N, param); % Anonymous handle to function

    reltol = 1E-6;
    options = odeset('RelTol', reltol);

    % Runge-Kutta implementation
    [tout, Nout] = ode45(odefun, tspan, N0, options);    

end

