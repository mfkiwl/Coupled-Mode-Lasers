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
M = 31;
deta = (etaAB_max - etaAB_min)/(M-1);
eta_AB = etaAB_min:deta:etaAB_max;

theta_max = 2*pi;
theta_min = 0.0;
N = 31;
dtheta = (theta_max - theta_min)/(N-1);
theta_n = theta_min:dtheta:theta_max;

[X, Y] = meshgrid(eta_AB, theta_n);

z = 0;
eig_vals = z*X.*Y;

etaAB = etaAB_min;
etaBA = etaAB_min;

[Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);

opt = 0;

for m = 1:M
    for n = 1:N
        
        etaAB = eta_AB(m);
        etaBA = eta_AB(m);
        theta = theta_n(n);
        
        [Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);
        
        eig_vals(n,m) = max(real(E));
        
    end
    
    disp(['Row: ' num2str(m)])
end

% Create colour map
values = eig_vals;
ncols = 21;

maxv = max(values(:));
minv = min(values(:));

vspan = maxv-minv;
dv = vspan/(ncols - 1);

sample = minv:dv:maxv;
neg_mask = (sample <= 0);
pos_mask = (sample > 0);

c0 = 0.4;

red = (1-c0)*abs(sample.*pos_mask/maxv) + c0*pos_mask;
green = (1-c0)*abs(sample.*neg_mask/minv) + c0*neg_mask;
blue = zeros(size(sample));

mymap = transpose([red; green; blue]);


figure('Name', 'Maximum real eigenvalues v coupling coefficient')
contourf(X, Y, eig_vals)

colormap(mymap)
colorbar
title([num2str(M) ' by ' num2str(N) ' grid'])
xlabel('\eta')
ylabel('\theta (rad)')

