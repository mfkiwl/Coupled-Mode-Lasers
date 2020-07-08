function [dN] = asymPairS(N0, param)
%ASYMPAIRS Provides the rate equations for the asymmetric double cavity
%model (for use with nonlinear solver)
%%
% *ASYMPAIRS*
%
%%  Description
%
% Provides the rate equations for the asymmetric double cavity model
%
%%  Usage
%
%    [dN] = asymPair(N0, param);
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
%       yn = param.yn           1/(tau_N) - carrier recombination rate
%       kp = param.kp           1/(2*tau_p) - cavity loss rate
%       aH = param.aH           Linewidth enhancement factor
%       QA = param.QA           Normalised pumping rate in guide A    
%       QB = param.QB           Normalised pumping rate in guide A 
%       etaAB = param.etaAB     Amplitude of coupling coefficient AB 
%       etaBA = param.etaBA     Amplitude of coupling coefficient BA
%       theta = param.theta     Phase of coupling coefficient 
%       DW = param.DW           Detuning between the cavity resonances
%
%%  Returns
%
%    dN         An array of the derivatives of the variables in N0
%
%% Notes
% 
%   This function is just an interface to call asymPair without the
%   parameter t

    dN = asymPair(0, N0, param);

end

