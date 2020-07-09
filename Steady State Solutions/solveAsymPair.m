function [Ns, found, E, esign] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, varargin)
% SOLVEASYMPAIR Finds the steady state solutions for the asymmetric pair model 
%%
% *SOLVEASYMPAIR*
%
%%  Description
%
% Finds the steady state solutions for the asymmetric pair model 
%
%%  Usage
%
%   [Ns, found, E, esign, Iout] = solveAsymPair(QA, QB, etaAB, etaBA, theta, DW, param, varargin)
%
%%  Arguments
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
%                   opt(1)  if > 0, the routine reports to the command
%                           window
%
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
%	Ns          The steady state solutions (if found) in the form above
%
%	found       Boolean indicating solutions found
%
%    E        	Eigenvalues of Jacobian
%
%    esign   	-1 indicates stable solution, 1 indicates unstable solution
%
%%  Dependencies
%
%   This routine calls:
%
%       <singleSlab.html asymPairS> - the model rate equations
%
%%

    % Flags
    no_initial = true;
    report = false;
    found = false;
    
    tol = 1E-6;         % Tolerance used for near zero-valued eigenvalues 
    varnum = 5;         % Number of independent variables
    
    if (nargin > 7)
        
        if (varargin{1} > 0)
            
            report = true;
            
        end
        
    end
    
    % Check for intial values
    if (nargin > 8)

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

    
    % Set additional required parameters
    param.QA = QA;
    param.QB = QB;
    param.etaAB = etaAB;
    param.etaBA = etaBA;
    param.theta = theta;
    param.DW = DW;
      
    if (report)
        
        runstr = 'Asymmetric Pair (coupling coefficients set, not calculated)'; 
        
        % Extract waveguide parameters
        k0 = param.k0;
        n1 = param.n1;
        n2 = param.n2;
        
        % Set local variables to report
        yn = param.yn;
        kp = param.kp;
        aH = param.aH;
        theta = param.theta;

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
        disp(['kp = ' num2str(kp) ' - the cavity loss rate (ns^-1)']);
        disp(['aH = ' num2str(aH) ' - linewidth enhancement factor']); 
        disp(['yn = ' num2str(yn) ' - recombination rate (ns^-1)']);
        disp(' ');
        disp('Individual guide parameters:');
        disp(['QA = ' num2str(QA) ' - total pump power in guide (1)']);   
        disp(['QB = ' num2str(QB) ' - total pump power in guide (2)']);
        disp(' '); 
        
    end

    if (no_initial)
        
        % Set default initial conditions
        N0 = zeros(5,1);

        % Approximate steady state solutions
        MA = 1.0;   % Carrier concentration in guide A
        MB = 1.0;   % Carrier concentration in guide A

        YA = sqrt(QA/MA - 1);	% Optical amplitude in guide A
        YB = sqrt(QB/MB - 1);	% Optical amplitude in guide A

        phi = pi;           % Phase between optical fields 

        N0(1) = MA;         % Carrier concentration in guide A
        N0(2) = MB;         % Carrier concentration in guide B
        N0(3) = YA;     	% Amplitude in guide A
        N0(4) = YB;     	% Amplitude in guide B
        N0(5) = phi;        % Relative phase between fields in A and B
        
    end

    % Anonymous handle to function including passed parameters   
    fun = @(x) asymPairS(x, param);    

    % MaxFunctionEvaluations: '100*numberOfVariables' (default)

    %options = optimoptions('fsolve','Display','off','PlotFcn',@optimplotfirstorderopt);
    %options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations', 2400);
    options = optimoptions('fsolve','Display','off');
    %options = optimoptions('fsolve','Display','iter', 'FunctionTolerance', 1E-4);
    %options = optimoptions('fsolve','Display','iter', 'StepTolerance', 1E-5);
    %options = optimoptions('fsolve','Display','iter', 'OptimalityTolerance', 0.04);
    %options = optimoptions('fsolve','Display','off', 'FunctionTolerance', 1E-4);

    [Ns,~,xstat,output,J] = fsolve(fun, N0, options);

    if (xstat > 0)

        found = true;

        % Find eigenvectors and eigenvalues of Jacobian
        [E] = eig(J);

        Etest = E;
        Ezero = false;

        % Check for zero eigenvalues
        for n = 1:length(E)

            if ((E(n) == 0.0)||(abs(E(n)) < tol))

                % Make these negative
                Etest(n) = -1.0;
                Ezero = true;

            end

        end

        Er = real(Etest);
        Emax = max(Er);

        if (Emax < 0.0)

            esign = -1; % Stable solution

        else

            esign = 1;  % Unstable solution

        end


        if (report)
            
            disp(' ')
            disp('Solution found: ')
            disp(['MA = ' num2str(Ns(1))])
            disp(['MB = ' num2str(Ns(2))])
            disp(['YA = ' num2str(Ns(3))])
            disp(['YB = ' num2str(Ns(4))])
            disp(['phi = ' num2str(Ns(5))])

            disp(' ');
            disp('Eigenvalues of Jacobian: ');

            len = length(E);

            for n = 1:len

                disp(['E' num2str(n) ' = ' num2str(real(E(n))) ' + ' num2str(imag(E(n))) ' i']);

            end

            disp(' ');

            if (esign < 0.0)

                disp('Stable solution');

                if (Ezero) 

                    disp(' ');
                    disp('There were zero or near zero-valued eigenvalues');

                end

            else

                disp('Unstable solution');

            end

        end

    else

        % *** NO SOLUTION FOUND ***
        if (report)

            disp('No solution found');
            disp(' ');
            disp('Current options:');
            disp(options);
            disp(' ');

            switch (xstat)
                case 0
                    disp(['Number of iterations exceeded options.MaxIterations' ...
                    ' or number of function evaluations exceeded']);
                    disp('options.MaxFunctionEvaluations');
                case -1
                    disp('Output function or plot function stopped the algorithm.');
                case -2
                    disp('Equation not solved. The exit message can have more information.');
                case -3
                    disp(['Equation not solved. Trust region radius became too small'...
                    '(trust-region-dogleg algorithm).']);
            end

            disp(' ');
            disp(['Iterations = ' num2str(output.iterations)]);
            disp(['Function evaluations = ' num2str(output.funcCount)]);
            disp(['Algorithm = ' output.algorithm]);
            disp(['Measure of first-order optimality = ' num2str(output.firstorderopt)]);
            disp(output.message);
            disp(' ');

        end

        E = zeros(varnum,1);
        esign = 0.0;
        Ns = N0;

    end

end






