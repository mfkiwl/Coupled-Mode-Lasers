function [dN] = coupled1D(~, N0, param)
%COUPLED1D Rate equations for the 1D coupled mode model (twin cavity)
%%
% *COUPLED1D*
%
%%  Description
%
% Provides the rate equations for the coupled mode model (twin cavity). The
% implementation is for the model given in [1]. This routine is used 
% directly by the Runge-Kutta solver (ode45) and indirectly (via
% coupled1DS) in conjunction with the nonlinear solver (fsolve).
%  
% [1] M.J. Adams et al, Phys. Rev. A 95(5), 053869 (2017)
%
%%  Usage
%
%    [dN] = coupled1D(t, N0, param);
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
    eta = param.eta;        % Amplitude of coupling coefficient 
    theta = param.theta;    % Phase of coupling coefficient 
    DW = param.DW;          % Detuning between the cavity resonances

    dN = zeros(size(N0));

    % dMA/dt:
    dN(1) = yn*(QA - MA*(1.0 + IA));

    % dMB/dt:
    dN(2) = yn*(QB - MB*(1.0 + IB));

    % dYA/dt:
    dN(3) = kp*(MA - 1.0)*YA - eta*YB*sin(theta + phi);

    % dYB/dt
    dN(4) = kp*(MB - 1.0)*YB - eta*YA*sin(theta - phi);

    % dphi/dt WARNING - THIS TERM APPEARS TO INDUCE NUMERICAL ERRORS (see
    % notes)
    dN(5) = aH*kp*(MA - MB) - DW + eta*((YA/YB)*cos(theta - phi) - (YB/YA)*cos(theta + phi));

end

