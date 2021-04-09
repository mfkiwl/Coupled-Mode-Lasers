function [dN] = asymPair(~, N0, param)
%ASYMPAIR Rate equations for the asymmetric double cavity model
%%
% *ASYMPAIR*
%
%%  Description
%
% Provides the rate equations for the asymmetric double cavity model
%
%%  Usage
%
%    [dN] = asymPair(t, N0, param);
%
%%  Arguments
%
%    t          an array of time steps (not used in the routine, but 
%               required for Runge-Kutta implementation).
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
%%  Notes
% When comparing with the single guide case by putting eta = 0, it is
% found that discrepancies arise even though the equations above are
% decoupled and should be equivalent to the single guide case. These
% discrepancies are in the amplitude and frequency of the relaxation
% oscillations. However, cutting off the end of Eq. 5 solves the problem,
% even though this last term is numerically zero. 
%   Moreover, when passing the time parameter t explictly to the routine,
% is found that the time intervals with the full expression for Eq. 5 are
% different to those in the single guide routine. This suggests that the
% ode45 routine that implements the Runge Kutte routine is dynamically
% changing the step size - presumably to compensate for numerical
% instability. This may be arising in the calculation of the last term of
% Eq. 5 where the optical amplitudes may be close to zero, even though this
% term is zero when eta = 0.
%
%%

    MA = N0(1);         % Carrier concentration in guide A
    MB = N0(2);         % Carrier concentration in guide A
    YA = N0(3);         % Amplitude in guide A
    YB = N0(4);         % Amplitude in guide B
    phi = N0(5);        % Relative phase between fields in A and B

    % Intensities (all values real)
    IA = YA*YA;
    IB = YB*YB;

    yn = param.yn;          % 1/(tau_N) - carrier recombination rate
    kp = param.kp;          % 1/(2*tau_p) - cavity loss rate
    aH = param.aH;          % Linewidth enhancement factor
    QA = param.QA;          % Normalised pumping rate in guide A    
    QB = param.QB;          % Normalised pumping rate in guide A 
    etaAB = param.etaAB;    % Amplitude of coupling coefficient AB
    etaBA = param.etaBA;    % Amplitude of coupling coefficient BA
    theta = param.theta;    % Phase of coupling coefficient 
    DW = param.DW;          % Detuning between the cavity resonances

    dN = zeros(size(N0));

    % dMA/dt:
    dN(1) = yn*(QA - MA*(1.0 + IA));

    % dMB/dt:
    dN(2) = yn*(QB - MB*(1.0 + IB));

    % dYA/dt:
    dN(3) = kp*(MA - 1.0)*YA - etaAB*YB*sin(theta + phi);

    % dYB/dt:
    dN(4) = kp*(MB - 1.0)*YB - etaBA*YA*sin(theta - phi);

    % dphi/dt: 
    dN(5) = aH*kp*(MA - MB) - DW + etaAB*(YA/YB)*cos(theta - phi)...
        - etaBA*(YB/YA)*cos(theta + phi);

end

