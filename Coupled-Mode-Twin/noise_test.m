%% Noise test
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

tsim1 = 0.4; %50;         % Simulation time (in units of 1/yn)
QA = 10;            % Normalised pump power in laser A
QB = 10;            % Normalised pump power in laser B
etaAB = 0.53383;    % Amplitude of coupling coefficient of B laser in dYA/dt
etaBA = 0.53383;    % Amplitude of coupling coefficient of A laser in dYB/dt
theta = 0.0;        % Phase of complex coupling
DW = 0;             % Frequency detuning

% opt = 0;            % Do NOT plot graphs, do NOT zero out noise 
opt = 1;            % DO plot graphs, do not zero out noise
% opt = 2;            % Do NOT plot graphs, DO zero out noise
%opt = 3;            % DO plot graphs, DO zero out noise            
%%
% Call the routine to calculate temporal dynamics. 

[t, N] = gaussianNoise(tsim1, QA, QB, etaAB, etaBA, theta, DW, param, opt);



