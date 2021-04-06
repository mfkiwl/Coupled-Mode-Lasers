function [dN] = coupled1DS(N0, param)
%COUPLED1DS Rate equations for coupled mode model (no time dependence).
%%
% *COUPLED1DS*
%
%%  Description
%
% Provides the rate equations for the coupled mode model (twin cavity), 
% implementation is for the model given in [1]. This routine is called
% without any time dependence for use with the nonlinear solver (fsolve).
%
% The 'S' at the end of the function name indicates that it is used
% primarily for finding the steady state solutions.
%  
% [1] M.J. Adams et al, Phys. Rev. A 95(5), 053869 (2017)
%
%%  Usage
%
%    [dN] = coupled1DS(N0, param);
%
%%  Arguments
%
%    N0         vector of inital conditions for the equations:
%
%    Carrier concentrations:
%       MA = N0(1)          Carrier concentration in guide A
%       MB = N0(2)          Carrier concentration in guide B
%    Optical fields:
%       YA = N0(3)          Amplitude in guide A
%       YB = N0(4)          Amplitude in guide B
%       phi = N0(5)         Relative phase between fields in A and B
%
%    param      an array containing the parameters:
%
%       yn = param.yn       1/(tau_N) - carrier recombination rate
%       kp = param.kp       1/(2*tau_p) - cavity loss rate
%       aH = param.aH       Linewidth enhancement factor
%       QA = param.QA       Normalised pumping rate in guide A    
%       QB = param.QB       Normalised pumping rate in guide A 
%       eta = param.eta     Amplitude of coupling coefficient 
%       theta = param.theta	Phase of coupling coefficient 
%       DW = param.DW       Detuning between the cavity resonances
%
%%  Returns
%
%    dN         An array of the derivatives of the variables in N0
%
%%  Dependencies
%
%   This routine calls:
%
%       coupled1D(t, N, param)      - the model rate equations
%
%%


    % Call time dependent routine with dummy variable
    dN = coupled1D(0, N0, param);

end

