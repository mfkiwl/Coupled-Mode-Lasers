%% Pulse simulation
%
% Simple implementation of pulse behaviour
%
%% Code
%
% First load the laser parameters (these include the width of the guide)

% param = loadParams('PRA_95_053869');  % Laser cavity parameters using in Ref [1]

% Step size (set to zero or less to use built-in routine)
DT = 0.001;

%% 
% The param structure contains the following fields. Uncomment and edit the
% line if you want to change any of the values
%
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

tsim1 = 50; % Simulation time (in units of 1/yn)
QA = 10;    % Normalised pump power in laser A
QB = 10;    % Normalised pump power in laser B
d = 16;     % Distance between guides (microns)
DW = 0;     % Frequency detuning

opt = 0;    % Do not plot graphs 
%%
% Call the routine to calculate temporal dynamics. Note that this routine
% calls the routine 'singleSlab', which calculates the coupling coefficent
% (real only) for the guides with the refractive index steps (n1 and n2), 
% width (w) andedge-to-edge separation (d) specified above.

[tout1, Nout1] = runCoupled1D_dt(tsim1, QA, QB, d, DW, param, opt, 0, DT); 

%%
% This runs the simulation for tsim*tau_N nanoseconds.
%
% Next, save the dynamic variables at the end of the simulation to use as
% initial values for the next run.

N0 = Nout1(end,:);   % End values of simulation
N0 = transpose(N0); % Transpose array to be in 5 x 1 form

%%
% Enter new simulation and pump values (edit these)

t_pulse = 3;           % Simulation time (short pulse)
QA_pulse = 20;      % Normalised pump power in laser A
QB_pulse = 5;       % Normalised pump power in laser B

%%
% Call runCoupled1D_dt again, this time with intial values and pulse pump

[tout2, Nout2] = runCoupled1D_dt(t_pulse, QA_pulse, QB_pulse, d, DW, param, opt, N0, DT);

%%
% Again, get end values for initial values of next run

N0 = Nout2(end,:);   % End values of simulation
N0 = transpose(N0); % Transpose array to be in 5 x 1 form

%%
% Enter new simulation time

tsim2 = tsim1 - t_pulse;    % Total simulation time will be 2*tsim1

% Call runCoupled1D_dt one more time with original pump values

[tout3, Nout3] = runCoupled1D_dt(tsim2, QA, QB, d, DW, param, opt, N0, DT);

%%
% Concatonate time and dynamic variables together

t1 = tout1(end);        % End time of tout1
tout2 = tout2 + t1;     % Add to tout2
t2 = tout2(end);        % End time of tout2
tout3 = tout3 + t2;     % Add to tout3

t = cat(1, tout1, tout2, tout3);
N = cat(1, Nout1, Nout2, Nout3);

%%
% Output values (uncomment any others you want to plot)
% MA  = N(:,1);
% MB  = N(:,2);
YA  = N(:,3);
YB  = N(:,4);
% phi = N(:,5);

%%
% Calculate intensities from optical fields
IA = YA.*YA;
IB = YB.*YB;

%%
% Plot results

lw = 1.5;       % graph linewidth
figure;
hold on;
plot(t, IA, 'LineWidth', lw);
plot(t, IB, 'LineWidth', lw);
title('Intensities');
ylabel('Optical intensity');
xlabel('Simulation time (\tau_{N})');
legend('I_{1}','I_{2}');
grid on;

%%
% Zoom in on pulse

xlim([47 56]);
