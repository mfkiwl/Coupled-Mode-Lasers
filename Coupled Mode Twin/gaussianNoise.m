%% Asymmetric laser pair model with Gaussian noise
function [tout, Nout] = gaussianNoise(tsim, QA, QB, etaAB, etaBA, theta, DW, param, varargin)
%
%% Description
%
%%  Usage
%
%   [tout, Nout] = runCoupled1D(tsim, QA, QB, d, DW, param, varargin)
%
%%  Arguments
%
%   tsim        integer specifying the simulation time tspan via tspan =
%               tsim/yn, where yn is the recombination rate (in 1/ns) 
%
%   QA          Pump power in guide (1)
%
%   QB          Pump power in guide (2)
%
%   etaAB       Amplitude of coupling coefficient of B laser in dYA/dt
%
%   etaBA       Amplitude of coupling coefficient of A laser in dYB/dt
%
%   theta       Phase of complex coupling 
%
%   DW          frequency detuning (rad/ns)
%
%   param       structure containing laser parameters:
%
%               param.aH    Linewidth enhancement factor
%               param.k0    free space wavevector (1/micron)
%               param.kp    1/(2*tau_p) - cavity loss rate (1/ns)
%               param.n1    core refractive index 
%               param.n2    cladding refractive index
%               param.yn    1/(tau_N) - carrier recombination rate (1/ns)
%
%   varargin    optional values:
%
%                   opt(1)  a numeric value. If greater than zero, the
%                           graphs are plotted
%                   opt(2)  a 5 x 1 array containing initial values:
%
%                       N0(1)   Carrier concentration in guide A
%                       N0(2)   Carrier concentration in guide B
%                       N0(3)   Optical amplitude in guide A
%                       N0(4)   Optical amplitude in guide B
%                       N0(5)   Relative phase between fields in A and B
%
%%  Returns
% 
%   Nout        4001 x 5 array containing the time evolution of variables
%
%                   Nout(:,1)   Carrier concentration in guide A
%                   Nout(:,2)   Carrier concentration in guide B
%                   Nout(:,3)   Optical amplitude in guide A
%                   Nout(:,4)   Optical amplitude in guide B
%                   Nout(:,5)   Relative phase between fields in A and B
%
%   tout        4001 x 1 array of corresponding time values
%
%% Notes
%
% List of all parameters:
%
%       k0          free space wavevector (1/micron)
%       n1          core refractive index 
%       n2          cladding refractive index
%       yn          1/(tau_N) - carrier recombination rate
%       kp          1/(2*tau_p) - cavity loss rate
%       aH          Linewidth enhancement factor
%       QA          Normalised pumping rate in guide A    
%       QB          Normalised pumping rate in guide A 
%       etaAB       Amplitude of coupling coefficient AB
%       etaBA       Amplitude of coupling coefficient BA
%       theta       Phase of coupling coefficient 
%       DW          Detuning between the cavity resonances
%
%% Code
%
% Code begins here

%% Check input arguments

    % Flag to plot graphs
    plotgraphs = false;
    zero = false;
    no_initial = true;

    % Check for flag to plot graphs
    if (nargin > 8)
        
        opt = varargin{1};
        
        switch opt
            
            case 0
                plotgraphs = false;
                zero = false;
            case 1
                plotgraphs = true;
                zero = false;
            case 2
                plotgraphs = false;
                zero = true;
            case 3
                plotgraphs = true;
                zero = true;
        end

    end
    
    % Check for intial values
    if (nargin > 9)

        N0 = varargin{2};
        
        sz = size(N0);
        
        if ((sz(1) ~= 5)||(sz(2)~=1))
            
            error('Initial values are not given in a 5 x 1 array')
            
        end
        
        no_initial = false;

    end

    if (QA < 0)
        error('eta1 cannot be negative')
    end

    if (QB < 0)
        error('eta2 cannot be negative')
    end

 
%% Set local variables
 
    % Extract waveguide parameters
    k0 = param.k0;
    n1 = param.n1;
    n2 = param.n2;
    
    % Set local variables
    yn = param.yn;
    kp = param.kp;
    aH = param.aH;
    
    
%% Set up Gaussian noise sources

    % Photon lifetime (ns) 
    tau_p = 1/(2.0*kp);

    % Multiple of photon lifetimes for the period of the fluctuations.
    Np = 32;

    % Period of the fluctuations
    DT = Np*tau_p;

    % Number of points in time arrays
    Npt = ceil(tsim/DT);

    if not(zero)
        
        % Create time arrays
        seed1 = 1;
        seed2 = 2;
        rng(seed1);
        xiA = randn(Npt, 1);
        rng(seed2);
        xiB = randn(Npt, 1);

        % Numerical factors
        beta_sp = 2.0e-5;       % Spontaneous emission factor
        c = 29.9792458;         % Speed of light (cm/ns)
        a_diff = 1e-15;         % Differential gain (cm^2)
        N_0 = 1e18;             % Transparency carrier concentration (cm^-3)
        ng = 3.4;               % Group refractive index
        Gam = 1;                % Optical confinement factor

        F = sqrt((beta_sp*c*a_diff/n1)*(N_0 + ng/(Gam*c*a_diff*tau_p))/DT);

        xiA = F*xiA;
        xiB = F*xiB;
    
    else
        
        xiA = zeros(Npt, 1); 
        xiB = zeros(Npt, 1);
        
    end
        
    xiA2 = sum(xiA.*xiA)/Npt;
    xiB2 = sum(xiB.*xiB)/Npt;
    xiAB = sum(xiA.*xiB)/Npt;    
    

%% Initial conditions

    if (no_initial)
        % Set default initial conditions
        N0 = zeros(5,1);
        %
        %   Carrier concentrations:
        %       MA = N0(1);         % Carrier concentration in guide A
        %       MB = N0(2);         % Carrier concentration in guide A
        %   Optical fields:
        %       YA = N0(3);         % Amplitude in guide A
        %       YB = N0(4);         % Amplitude in guide B
        %       phi = N0(5);        % Relative phase between fields in A and B

        MA_in = 0.0;   % Guide 1
        MB_in = 0.0;   % Guide 2

        % Steady state solutions
        % dM = 0.000001;
        % MS = 1.0;
        % MA = MS + dM;   % Guide 1
        % MB = MS + dM;   % Guide 2

        IA_in = 10;%1E-4;           % Initial signal power symmetric mode 
        IB_in = 10;%1E-4;           % Initial signal power anti-symmetric mode 
        YA_in = sqrt(IA_in);    % optical amplitude in A
        YB_in = sqrt(IB_in);	% optical amplitude in A 

        % Steady state solutions
        % dA = 0.000001;
        % YAS = real(sqrt(QA - 1));
        % YBS = real(sqrt(QB - 1));
        % YA = YAS + dA;   % Guide 1
        % YB = YBS + dA;   % Guide 2

        phi_in = pi; %0.001;        % Phase 

        N0(1) = MA_in;	% Carrier concentration in guide A
        N0(2) = MB_in;	% Carrier concentration in guide B
        N0(3) = YA_in;	% Amplitude in guide A
        N0(4) = YB_in;	% Amplitude in guide B
        N0(5) = phi_in;	% Relative phase between fields in A and B
        
    end
    
%% Set up and run Runge-Kutta routine

    % Time span (yn = 1/tau, where tau is the lifetime)
    maxt = tsim/yn;
    mint = 0.0;
    npts = 4001; 
    dt = (maxt - mint)/(npts - 1.0);

    tspan = mint:dt:maxt;

    odefun = @(t, N) asymPairGauss(t, N); % Anonymous handle to function

    reltol = 1E-6;
    options = odeset('RelTol', reltol);
    
    % Runge-Kutta implementation
    lastwarn('');
    
    [tout, Nout] = ode45(odefun, tspan, N0, options);  
    
    [warnMsg, ~] = lastwarn;
    
    % TEST
    tpts = length(tout);
    dNout = zeros(size(Nout));
    
    for n = 1:tpts
        
        if (n == tpts)
            
            disp(tout(n))
            
        end
        
        dNout(n,:) = asymPairGauss(tout(n), Nout(n,:));
         
    end
        
        

 
%% Report input parameters

    runstr = 'Asymmetric Pair with Gaussian Noise (coupling coefficients set, not calculated)'; 

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
    disp(['etaAB = ' num2str(etaAB) ' - Amplitude of coupling coefficient of B laser in dYA/dt']);
    disp(['etaBA = ' num2str(etaBA) ' - Amplitude of coupling coefficient of A laser in dYB/dt']);
    disp(['theta = ' num2str(theta) ' - coupling coefficient (argument)']);
    disp(['DW = ' num2str(DW) ' - frequency detuning']);

    disp(' ');
    disp('General parameters:');
    disp(['tau_p = ' num2str(tau_p) ' - photon lifetime (ns)'])
    disp(['kp = ' num2str(kp) ' - the cavity loss rate (ns^-1)']);
    disp(['aH = ' num2str(aH) ' - linewidth enhancement factor']); 
    disp(['yn = ' num2str(yn) ' - recombination rate (ns^-1)']);
    disp(' ');
    disp('Individual guide parameters:');
    disp(['QA = ' num2str(QA) ' - total pump power in guide (1)']);   
    disp(['QB = ' num2str(QB) ' - total pump power in guide (2)']);
    disp(' '); 
    disp(['Simulation time = ' num2str(tsim) '/yn = ' num2str(tsim/yn) ' ns']);
    disp(' ');
    disp('Gaussian noise generation:')
    disp(['Number of noise terms: ' num2str(Npt)]);
    disp(['variance of (normalised) noise in A: ' num2str(xiA2)]);
    disp(['variance of (normalised) noise in B: ' num2str(xiB2)]);
    disp(['covariance of (normalised) noise between A and B: ' num2str(xiAB)]);
    
    
%% Report results

    % Output values
    MAt = Nout(:,1);
    MBt = Nout(:,2);
    YAt = Nout(:,3);
    YBt = Nout(:,4);
    phit = Nout(:,5);

    % Test derivatives at end of calculation
    Nt = transpose(Nout(end,:));
    dN = asymPairGauss(maxt, Nt);

    dMAdt = dN(1);
    dMBdt = dN(2);
    dYAdt = dN(3);
    dYBdt = dN(4);
    dphidt = dN(5);

    % Output results
    disp(' ');
    disp('At end of simulation:');
    disp(['N(1) = ' num2str(MAt(end))]);
    disp(['N(2) = ' num2str(MBt(end))]);
    disp(['A_{1} = ' num2str(YAt(end))]);
    disp(['A_{2} = ' num2str(YBt(end))]);
    disp(['phi = ' num2str(phit(end))]);
    disp(' ');

    disp('Derivatives:');
    disp(['dN1/dt = ' num2str(dMAdt)]);
    disp(['dN2/dt = ' num2str(dMBdt)]);
    disp(['dA1/dt = ' num2str(dYAdt)]);
    disp(['dA2/dt = ' num2str(dYBdt)]);
    disp(['dphi/dt = ' num2str(dphidt)]);
    
    if not(isempty(warnMsg))
        
        disp('')
        disp('')
        disp('THE ROUTINE ENCOUNTERED A WARNING:')
        disp('')
        disp(warnMsg)
        disp('')
        disp(['Time at end of simulation: ' num2str(tout(end))])
        
    end

%% Optionally plot graphs

	if plotgraphs
        
        % Plot noise
        test_t = 0:(tsim/(Npt-1)):tsim;
    
%         figure
%         plot(test_t, xiA)
%         grid on
%         title(['\xi_{A} (<\xi_{A}|\xi_{A}> = ' num2str(xiA2) ')'])
%         xlabel('time (ns)')
% 
%         figure
%         plot(test_t, xiB)
%         grid on
%         title(['\xi_{B} (<\xi_{B}|\xi_{B}> = ' num2str(xiB2) ')'])
%         xlabel('time (ns)')
 
        figure
        plot(test_t, xiA, test_t, xiB)
        grid on
        title(['\Delta T = ' num2str(DT) ' ns'])
        legend('\xi_{A}', '\xi_{B}')
        xlabel('time (ns)')

        % Plot solution
        lw = 1.5; % linewidth

        % Plot carriers 
        figure;
        hold on;
        plot(tout, MAt, 'LineWidth', lw);
        plot(tout, MBt, 'LineWidth', lw);
        title(['Carriers: QA = ' num2str(QA) ', QB = ' num2str(QB) ' (coupled)']);
        ylabel('Carrier concentration');
        xlabel('Simulation time (1/\gamma)');
        legend('N^{(1)}', 'N^{(2)}');
        grid on;
        
        % Plot amplitudes 
        figure;
        hold on;
        plot(tout, YAt, 'LineWidth', lw);
        plot(tout, YBt, 'LineWidth', lw);
        title(['Amplitudes: QA = ' num2str(QA) ', QB = ' num2str(QB) ' (coupled)']);
        ylabel('Amplitude');
        xlabel('Simulation time (1/\gamma)');
        legend('Y^{(1)}', 'Y^{(2)}');
        grid on;
        
        % Plot derivatives of amplitudes 
        figure;
        hold on;
        plot(tout, dNout(:,3), 'LineWidth', lw);
        plot(tout, dNout(:,4), 'LineWidth', lw);
        title(['dY/dt: QA = ' num2str(QA) ', QB = ' num2str(QB) ' (coupled)']);
        ylabel('dY/dt');
        xlabel('Simulation time (1/\gamma)');
        legend('dY^{(1)}/dt', 'dY^{(2)}/dt');
        grid on;
        
        IAt = YAt.*YAt;
        IBt = YBt.*YBt;

        % Plot intensities
        figure;
        hold on;
        plot(tout, IAt, 'LineWidth', lw);
        plot(tout, IBt, 'LineWidth', lw);
        title(['Intensities: QA = ' num2str(QA) ', QB = ' num2str(QB) ' (coupled)']);

        ylabel('Optical intensity');
        xlabel('Simulation time (1/\gamma)');
        legend('I_{1}','I_{2}');
        grid on;

        % Plot phase
        figure;
        hold on;
        plot(tout, phit, 'LineWidth', lw);
        title(['Phase: eta1 = ' num2str(QA) ', eta2 = ' num2str(QB)]);
        ylabel('Relative phase');
        xlabel('Simulation time (1/\gamma)');
        %ylim([-2*pi 2*pi]);
        grid on;

	end


%% Local rate equation function

    function [dN] = asymPairGauss(t, N0) 
    
    %%  Description
    %
    % Provides the rate equations for the asymmetric double cavity model
    %
    %%  Usage
    %
    %    [dN] = asymPairGauss(t, N0);
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
        
        % Noise terms
        index = floor(t/DT) + 1;
        
%         GA = xiA(index);
%         GB = xiB(index);
        
        % TEST
        GA = abs(xiA(index));
        GB = abs(xiB(index));
        
        dN = zeros(size(N0));

        % dMA/dt:
        dN(1) = yn*(QA - MA*(1.0 + IA));

        % dMB/dt:
        dN(2) = yn*(QB - MB*(1.0 + IB));

        % dYA/dt:
        dN(3) = kp*(MA - 1.0)*YA - etaAB*YB*sin(theta + phi) + GA;

        % dYB/dt:
        dN(4) = kp*(MB - 1.0)*YB - etaBA*YA*sin(theta - phi) + GB;

        % dphi/dt: 
        % TEST
        %dN(5) = 0;
        dN(5) = aH*kp*(MA - MB) - DW + etaAB*(YA/YB)*cos(theta - phi)...
            - etaBA*(YB/YA)*cos(theta + phi);

    end

end