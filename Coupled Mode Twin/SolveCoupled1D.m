function [Ns, found, E, esign] = SolveCoupled1D(QA, QB, d, DW, param, varargin)
% SolveCoupled1D Finds steady state solutions for coupled mode model
%
% Usage:
%
%   [Ns, found, E, esign] = SolveCoupled1D(QA, QB, d, DW, param, varargin);
%
% Arguments:
%
%   QA          Pump power in guide (1)
%
%   QB          Pump power in guide (2)
%
%   d           distance between guides (um)
%
%   DW          frequency detuning
%
%   param       structure containing parameters
%
% Return values:
%
%   esign       Sign of largest eigenvalue of Jacobian if -1 steady state
%               solution is stable, if 1 it is unstable
%
%   Ns          Steady state solutions
%
%   E           Eigenvalues of Jacobian
%
%   found       True if solution found, false otherwise
%

found = false;
report = false;
esign = 0;

if (QA < 0)
    error('eta1 cannot be negative')
end

if (QB < 0)
    error('eta2 cannot be negative')
end

if (nargin > 5)
    
    report = true;
    
end

% Extract waveguide parameters
k0 = param.k0;
n1 = param.n1;
n2 = param.n2;
w = param.w;

% Find coupling coefficient
[~, eta, ~, neff] = singleInt(k0, n1, n2, w, d);

% Set additional required parameters
param.eta = abs(eta);
param.QA = QA;
param.QB = QB;
param.DW = DW;
param.theta = phase(eta);

% Set local variables to report
yn = param.yn;
kp = param.kp;
aH = param.aH;
theta = param.theta;

if (report)

    runstr = 'Coupled cavity model'; 

    % Output parameters 
    disp(' ');
    disp(datestr(now));
    disp(' ');
    disp(runstr);
    disp(' ');
    disp('Guide parameters:');
    disp(['k0 = ' num2str(k0) ' - free-space wavevector (rad/micron)']);
    disp(['n1 = ' num2str(n1) ' - refractive index in guides']);
    disp(['n2 = ' num2str(n2) ' - refractive index outside guides']);
    disp(['neff = ' num2str(neff) ' - effective index of mode']);
    disp(['d = ' num2str(d) ' - distance between guides (microns)']);
    disp(['w = ' num2str(w) ' - width of guides (microns)']);
    disp(['|eta| = ' num2str(eta) ' - coupling coefficient (modulus)']);
    disp(['theta = ' num2str(theta) ' - coupling coefficient (argument)']);
    disp(['DW = ' num2str(DW) ' - frequency detuning']);

    disp(' ');
    disp('General parameters:');
    disp(['kp = ' num2str(kp) ' - the cavity loss rate (ns^-1)']);
    disp(['aH = ' num2str(aH) ' - linewidth enhancement factor']); 
    disp(['yn = ' num2str(yn) ' - recombination rate (ns^-1)']);
    disp(' ');
    disp('Individual guide parameters:');
    disp(['QA = ' num2str(QA) ' - total pump power in guide (1)']);   
    disp(['QB = ' num2str(QB) ' - total pump power in guide (2)']);
    

end


% Set initial conditions
x0 = zeros(5,1);
%
%   Carrier concentrations:
%       x0(1) = MA;         % Carrier concentration in guide A
%       x0(2) = MB;         % Carrier concentration in guide A
%   Optical fields:
%       x0(3) = YA;         % Amplitude in guide A
%       x0(4) = YB;         % Amplitude in guide B
%       x0(5) = phi;        % Relative phase between optical fields


% Approximate steady state solutions (for isolated guide)
dM = 0.1;
MS = 1.0;
MA = MS + dM;   % Guide 1
MB = MS + dM;   % Guide 2

% Approximate steady state solutions (for isolated guide)
dA = 0.001;
YAS = real(sqrt(QA - 1));
YBS = real(sqrt(QB - 1));
YA = YAS + dA;   % Guide 1
YB = YBS + dA;   % Guide 2

% Create phase shift
dphi = 0.01;
phiS = pi; 
phi = phiS + dphi;

if (report)
    disp(' ');
    disp('Predicted steady state solution: ');
    disp(['MAS = ' num2str(MS)]);
    disp(['MBS = ' num2str(MS)]);
    disp(['|YAS| = ' num2str(YAS)]);
    disp(['|YBS| = ' num2str(YBS)]);
    disp(['IAS = ' num2str(YAS*YAS)]);
    disp(['IAS = ' num2str(YBS*YBS)]);
    disp(['phi = ' num2str(phiS)]);
end

x0(1) = MA;         % Carrier concentration in guide A
x0(2) = MB;         % Carrier concentration in guide B
x0(3) = YA;     	% Amplitude in guide A
x0(4) = YB;     	% Amplitude in guide B
x0(5) = phi;        % Relative phase between optical fields

% Anonymous handle to function   
fun = @(x) coupled1DS(x, param);    

options = optimoptions('fsolve','Display','off');

[Ns,fval,xstat,out,J] = fsolve(fun, x0, options);

if (xstat > 0)
    
    found = true;
    
    MAS = Ns(1);
    MBS = Ns(2);
    YAS = Ns(3);
    YBS = Ns(4);
    phi = Ns(5);
    
    if (report)
    
        disp(' ');
        disp('Found steady state solution: ');
        disp(['MAS = ' num2str(MAS)]);
        disp(['MBS = ' num2str(MBS)]);
        disp(['|YAS| = ' num2str(abs(YAS))]);
        disp(['|YBS| = ' num2str(abs(YBS))]);
        disp(['phi = ' num2str(phi)]);
        
    end
    
    % Find eigenvectors and eigenvalues of Jacobian
    [E] = eig(J);
    
    Er = real(E);
    Emax = max(Er);
    if (Emax < 0.0)
        
        esign = -1; % Stable solution
        
    else
        
        esign = 1;  % Unstable solution
        
    end

    % Eigenvectors in diagonal matrix
    E1 = E(1);
    E2 = E(2);
    E3 = E(3);
    E4 = E(4);
    E5 = E(5);
    
    if (report)
    
        disp(' ');
        disp('Eigenvalues of Jacobian: ');
        disp(['E1 = ' num2str(E1)]);
        disp(['E2 = ' num2str(E2)]);
        disp(['E3 = ' num2str(E3)]);
        disp(['E4 = ' num2str(E4)]);
        disp(['E4 = ' num2str(E5)]);
        
        if (esign < 0.0)
            
            disp('Stable solution');
            
        else
            
            disp('Unstable solution');
            
        end
    
    end
    
else
    
    if (report)
        
        disp('Fsolve did not find a solution.');
        disp(out.message);
        disp(['fval(1) = ' num2str(fval(1))]);
        disp(['fval(2) = ' num2str(fval(2))]);
        disp(['fval(3) = ' num2str(fval(3))]);
        disp(['fval(4) = ' num2str(fval(4))]);
        disp(['fval(5) = ' num2str(fval(5))]);
        
    end
    
    E = zeros(4,1);
    
end


end






