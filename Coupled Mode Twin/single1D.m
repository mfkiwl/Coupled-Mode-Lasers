function [dN] = single1D(~, N0, param)
%single1D Implementation of 1D single guide laser for solution via the 
%Runge Kutta model 
%   
%Synatax:
%
%       [dN] = single1D(t, N0, param);
%
%   t is an array of time steps (not used in the routine, but required for
%   Runge-Kutta implementation).
%
%   N0 is a vector of inital conditions for the equations:
%
%       N = N0(1);         % Carrier concentration 
%
%       A = N0(2);         % Optical amplitude 
%
%   param is an array containing the parameters:
%
%       yn = param.yn;      % 1/(tau_N) - carrier recombination rate
%       kp = param.kp;      % 1/(2*tau_p) - cavity loss rate
%       Q = param.Q;        % Normalised pumping rate in guide A    

    N = N0(1);      % Carrier concentration
    A = N0(2);      % Optical amplitude

    % Intensities 
    I = conj(A)*A;

    yn = param.yn;      % 1/(tau_N) - carrier recombination rate
    kp = param.kp;      % 1/(2*tau_p) - cavity loss rate
    Q = param.Q;        % Normalised pumping rate

    dN = zeros(size(N0));

    % dN/dt:
    dN(1) = yn*(Q - N*(1.0 + I));

    % dA/dt:
    dN(2) = kp*(N - 1.0)*A;


end

