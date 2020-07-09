%% Stability map Delta f v d/a
%
%% Parameters

param.aH = 2;       	% Linewidth enhancement factor
param.k0 = 4.8332;   	% Free-space wavevector (1/micron)
param.kp = 326.7974;	% Cavity loss rate 1/(2*tau_p) (1/ns)
param.n1 = 3.400971;	% Refractive index in core
param.n2 = 3.4;         % Refractive index in cladding
param.w = 8;            % Cavity width (micron)
param.yn = 1;           % 1/tau_N, where tau_N is the carrier lifetime

theta = 0;              % coupling coefficient real

%% Pump power
P = 2.0;                  % P/Pth 
CQ = 11.4;
Q = CQ*(P - 1) + P;     % Normalised pump power in both guides
QA = Q;
QB = Q;

%% Set up array for d/a
NG = 201; % 31

dmax = 2.5;
dmin = 1.0;
M = NG;
dd = (dmax - dmin)/(M - 1);

% Note - routine uses width of guide and edge to edge separation, so here
% we use a = param.w
d_over_a = dmin:dd:dmax;
eta = realEta(d_over_a, param.w, param.k0, param.n1, param.n2);


%% Set up array for Deltaf
maxF = 5.0;
minF = -5.0;
N = NG;
dF = (maxF - minF)/(N - 1);

DeltaF = minF:dF:maxF;

% Convert to angular frequency
dOmega = 2.0*pi*DeltaF;


%% Calculate map

[X, Y] = meshgrid(d_over_a, DeltaF);

z = 0;
eig_vals = z*X.*Y;

DW = dOmega(1);
etaAB = eta(1);
etaBA = eta(1);

opt = 0;

[N0, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);

opt = 0;

h = waitbar(0, 'Generating stability map...');

for n = 1:N 
    
    for m = 1:M
        
        if (m == 1)
            
            % Get solution from previous row
            Ns = N0;
            
        end
    
        etaAB = eta(m);
        etaBA = eta(m);
        DW = dOmega(n);
        
        [Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt, Ns);
        
        if (found)
        
            eig_vals(n,m) = max(real(E));
            
            if (n == 1)
                
                % Store for next row
                N0 = Ns;
                
            end
            
        else
            
            eig_vals(n,m) = NaN;
            
            for k = m:M
                
                eig_vals(n,k) = NaN;
                
            end
            
            break
                 
        end 
        
    end
    
    waitbar(n/N, h);
    
end

close(h)

%% Create colour map
values = eig_vals;
ncols = 11;

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

figure('Name', 'Delta f v d/a')
contourf(X, Y, eig_vals,'EdgeColor', 'none')
colormap(mymap)
colorbar
title(['P/P_{th} = ' num2str(P)])
xlabel('{\it d/a}')
ylabel('\Delta{\it f} (GHz)')
grid on

% Clean up
clear

