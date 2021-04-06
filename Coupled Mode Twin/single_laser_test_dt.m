%% Single laser simulation to test fixed step Runge Kutta routine

% This examines the effect of different step sizes for single laser transients

% Use this code in any implementations using fixed step RK routine (RK4)
% Calls RK4 - check if this is present on the path!
if (exist('RK4','file') ~= 2)

    % Folder containing RK4.m
    addpath([userpath '/Coupled Mode Lasers/Numeric']); 

end


%% Set parameters

% Simulation time (multiple of 1/yn = tau_N)
tsim = 1.5;

% Carrier decay rate (1/ns)
yn = 1.0;

% Cavity loss rate (1/ns)
kp = 326;

% Normalised pump power
Q = 10;

% Assign parameter structure
param.yn = yn;
param.kp = kp;
param.Q = Q;

%% Step sizes for fixed-step RK routine

DT1 = 0.01;     
DT1 = DT1/yn;   

DT2 = 0.008;     
DT2 = DT2/yn;  

DT3 = 0.004;     
DT3 = DT3/yn; 

DT4 = 0.001;     
DT4 = DT4/yn; 

% Initial conditions
N0 = zeros(2,1);
N0(1) = 0.0;    % Carrier concentration
N0(2) = 1E-4;   % Optical amplitude

%% Run simulation using MATLAB RK solver

% Time span (yn = 1/tau, where tau is the lifetime)
t1 = tsim/yn;
t0 = 0.0;

odefun = @(t, N) single1D(t, N, param); % Anonymous handle to function


%% Run simulation using fixed-step RK solver
disp(' ')
disp(['Running with DT = ' num2str(DT1)])
tic
[tout1, Nout1] = RK4(odefun, N0, t0, t1, DT1);
toc

disp(' ')
disp(['Running with DT = ' num2str(DT2)])
tic
[tout2, Nout2] = RK4(odefun, N0, t0, t1, DT2);
toc

disp(' ')
disp(['Running with DT = ' num2str(DT3)])
tic
[tout3, Nout3] = RK4(odefun, N0, t0, t1, DT3);
toc

disp(' ')
disp(['Running with DT = ' num2str(DT4)])
tic
[tout4, Nout4] = RK4(odefun, N0, t0, t1, DT4);
toc

lw = 1.5; % Graph line-width for plotting

figure
hold on
plot(tout1, Nout1(:,1), 'LineWidth', lw)
plot(tout2, Nout2(:,1), 'LineWidth', lw)
plot(tout3, Nout3(:,1), 'LineWidth', lw)
plot(tout4, Nout4(:,1), 'LineWidth', lw)
title('Carrier concentration')
legend(['DT = ' num2str(DT1)], ['DT = ' num2str(DT2)],...
    ['DT = ' num2str(DT3)], ['DT = ' num2str(DT4)])
grid on
ylim([-1 2.5])
xlim([0 1])

figure
hold on
plot(tout1, Nout1(:,2), 'LineWidth', lw)
plot(tout2, Nout2(:,2), 'LineWidth', lw)
plot(tout3, Nout3(:,2), 'LineWidth', lw)
plot(tout4, Nout4(:,2), 'LineWidth', lw)
title('Optical amplitude')
legend(['DT = ' num2str(DT1)], ['DT = ' num2str(DT2)],...
    ['DT = ' num2str(DT3)], ['DT = ' num2str(DT4)])
grid on

ylim([-2 18])
xlim([0 1])


    
    