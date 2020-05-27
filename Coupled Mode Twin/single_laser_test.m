%% Single laser simulation to test fixed step Runge Kutta routine

% Calls RK4 - check if this is present on the path!
if (exist('RK4','file') ~= 2)

    % Folder containing bisect.m
    addpath([userpath '/Overlap-Factor-Model/Numeric']); 

end

%% Test of fixed-step RK routine
DT = 0.001;      % Edit this...
DT = DT/yn;     % ... not this!

%% Run simulation using MATLAB RK solver

% Simulation time (multiple of 1/yn)
tsim = 1.5;

% Carrier decay rate (1/ns)
yn = 1.0;

% Cavity loss rate (1/ns)
kp = 326;

% Normalised pump power
Q = 10;

% Set parameters
param.yn = yn;
param.kp = kp;
param.Q = Q;

% Initial conditions
N0 = zeros(2,1);
N0(1) = 0.0;    % Carrier concentration
N0(2) = 1E-4;   % Optical amplitude

% Time span (yn = 1/tau, where tau is the lifetime)
t1 = tsim/yn;
t0 = 0.0;
npts = 4001; 
dt = (t1 - t0)/(npts - 1.0);

tspan = t0:dt:t1;

odefun = @(t, N) single1D(t, N, param); % Anonymous handle to function

options = odeset('MaxStep', DT);
%options = odeset('MaxStep', 20*DT);

% Runge-Kutta implementation
[tout1, Nout1] = ode45(odefun, tspan, N0, options); 




[tout2, Nout2] = RK4(odefun, N0, t0, t1, DT);

lw = 1.5; % Graph line-width for plotting

figure
hold on
plot(tout1, Nout1(:,1), 'LineWidth', lw)
plot(tout2, Nout2(:,1), 'LineWidth', lw)
title(['Carrier concentration (DT = ' num2str(DT) ' ns)'])
legend('Variable step', 'Fixed step')
grid on

figure
hold on
plot(tout1, Nout1(:,2), 'LineWidth', lw)
plot(tout2, Nout2(:,2), 'LineWidth', lw)
title(['Optical amplitude (DT = ' num2str(DT) ' ns)'])
legend('Variable step', 'Fixed step')
grid on





    
    