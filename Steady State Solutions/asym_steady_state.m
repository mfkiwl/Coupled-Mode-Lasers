%% Pulse simulation
%
% Simple implementation of pulse behaviour
%
%% Code
%
% First load the laser parameters (these include the width of the guide)

%% 
% The param structure contains the following fields. Edit the
% line if you want to change any of the values

param.aH = 2;             % Linewidth enhancement factor
param.k0 = 4.8332;        % Free-space wavevector (1/micron)
param.kp = 326.7974;      % Cavity loss rate 1/(2*tau_p) (1/ns)
param.n1 = 3.400971;      % Refractive index in core
param.n2 = 3.4;           % Refractive index in cladding
param.w = 8;              % Cavity width (micron)
param.yn = 1;             % 1/tau_N, where tau_N is the carrier lifetime

%%
% Next, set the operating conditions and simulation time (edit this as
% necessary)

QA = 10;            % Normalised pump power in laser A
QB = 10;            % Normalised pump power in laser B
% etaAB = 0.53383;    % Amplitude of coupling coefficient of B laser in dYA/dt
% etaBA = 0.53383;    % Amplitude of coupling coefficient of A laser in dYB/dt
theta = 0.0;        % Phase of complex coupling
DW = 0;             % Frequency detuning

opt = 1;            % Report to the command window
%%

% Stability map of DW v eta
etaAB_max = 4.0;
etaAB_min = 0.0;
M = 11;
deta = (etaAB_max - etaAB_min)/(M-1);
eta_AB = etaAB_min:deta:etaAB_max;

etaBA_max = 4.0;
etaBA_min = 0.0;
N = 11;
deta = (etaBA_max - etaBA_min)/(N-1);
eta_BA = etaAB_min:deta:etaAB_max;

[X, Y] = meshgrid(eta_AB, eta_BA);

z = 0;
Z = z*X.*Y;

etaAB = etaAB_min;
etaBA = etaBA_min;

[Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);

opt = 0;

for m = 1:M
    for n = 1:N
        
        etaAB = eta_AB(m);
        etaBA = eta_BA(n);
        
        disp(['(' num2str(m) ', ' num2str(n) ')'])
        
        [Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);
        
        Z(n,m) = esign;
        %Z(n,m) = max(real(E));
        
    end
end

contourf(X, Y, Z);

