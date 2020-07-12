%% Stability map Pth v d/a
%
%% Parameters

param.aH = 2;       	% Linewidth enhancement factor
param.k0 = 4.8332;   	% Free-space wavevector (1/micron)
param.kp = 326.7974;	% Cavity loss rate 1/(2*tau_p) (1/ns)
param.n1 = 3.400971;	% Refractive index in core
param.n2 = 3.4;         % Refractive index in cladding
param.w = 8;            % Cavity width (micron)
param.yn = 1;           % 1/tau_N, where tau_N is the carrier lifetime

DW = 0;                 % zero detuning
theta = 0;              % coupling coefficient real

%% Set up array for d/a
NG = 31;

dmax = 2.5;
dmin = 1.0;
M = NG;
dd = (dmax - dmin)/(M - 1);

% Note - routine uses width of guide and edge to edge separation, so here
% we use a = param.w
d_over_a = dmin:dd:dmax;
eta = realEta(d_over_a, param.w, param.k0, param.n1, param.n2);


%% Set up array for Pth
maxP = 3.5;
minP = 1.001;
N = NG;
dP = (maxP - minP)/(N - 1);

P = minP:dP:maxP;
CQ = 11.4;
Q = CQ*(P - 1) + P;

%% Calculate map

[X, Y] = meshgrid(d_over_a, P);

z = 0;
eig_vals = z*X.*Y;

QA = Q(1);
QB = Q(1);
etaAB = eta(1);
etaBA = eta(1);

opt = 0;

[N0, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt);

if not(found)
    
    disp(' ')
    disp('NO SOLUTION FOUND FOR:')
	disp(['d/a = ' num2str(d_over_a(1)) ', ' 'Q = ' num2str(Q(1))])
    
end

opt = 0;

h = waitbar(0, 'Generating stability map...');

for m = 1:M
    
    for n = 1:N
        
        if (n == 1)
            
            % Get solution from previous row
            Ns = N0;
            
        end
    
        etaAB = eta(m);
        etaBA = eta(m);
        QA = Q(n);
        QB = Q(n);
        
        [Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, opt, Ns);
        
        if (found)
        
            eig_vals(n,m) = max(real(E));
            
            if (n == 1)
                
                % Store for next row
                N0 = Ns;
                
            end
            
        else
            
            eig_vals(n,m) = NaN;
            
            disp(' ')
            disp('NO SOLUTION FOUND FOR:')
            disp(['d/a = ' num2str(d_over_a(m)) ', ' 'Q = ' num2str(Q(n))])
                 
        end 
        
    end
    
    waitbar(m/M, h);
    
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


figure('Name', 'P/Pth v d/a')
contourf(X, Y, eig_vals,'EdgeColor', 'none')
colormap(mymap)
colorbar
xlabel('{\it d/a}')
ylabel('{\it P/P_{th}}')

Z = eig_vals;

save('Pth.mat', 'X', 'Y', 'Z');

% Clean up
%clear

