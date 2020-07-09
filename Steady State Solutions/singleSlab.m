function [GAA, eta, mode, neff] = singleSlab(k0, n1, n2, w, d, varargin)
%SINGLESLAB Calculates parameters for a single slab guide in the coupled mode model. 
%%
% *SINGLESLAB*
%
%%  Description
%
% This routine calculates parameters for a single guide in the coupled mode 
% model. Only the lowest even solution found is considered. In particular,
% the routine calculates the optical confinement factor, the effective
% index and the coupling coefficient for equal guides separated by a
% distance d. 
%
% The routine also returns the spatial profile of the mode and optionally
% plot this.
% 
%%  Usage
%
%   [GAA, eta, mode, neff] = singleInt(k0, n1, n2, w, d, varargin)
%
%%  Arguments
%
%   k0          free space wavevector (1/micron)
%   n1          core refractive index
%   n2          cladding refractive index
%   w           guide width (micron)
%   d           separation between guides (micron)
%   varargin    optional argument. If passed and greater than 0, the
%               routine reports to the command line and plots a graph of
%               the mode.
%
%%  Returns
%
%   GAA         the optical confinement factor for the guide
%   eta         the coupling coefficient for equal guides separated by d
%   mode        2 x 10001 array containing:
%                   mode(1,:) - spatial coordinates
%                   mode(2,:) - mode profile
%   neff        the effective refractive index
%
%%  Dependencies
%
%   This routine calls:
%
%       bisect          bisection routine*
%       trapezoid1D     1D integration routine*
%   
%   *Note that both bisect and trapezoid are in the 'Numeric' directory.
%   The routine checks to see if this are on the path and if not, attempts
%   to add it. For this to work, the 'Numeric' directory MUST be in a
%   folder called 'Photonic-Neurons', which should be on the MATLAB path
%
%       'C:\Users\<user_name>\Documents\MATLAB'  
%
%   or similar. Type 'userpath' at the command line to check.
%
%%

    % Boolean to indicate whether routine reports and plots results
    report = false;
    xspan = 9.0*w;
    c = 299792458E-3; % micron/ns

    precision = 1E-9;

    if (nargin > 5)

        if (varargin{1} > 0)
            report = true;
        end


    end

    % If profile is anti-guiding, report an error
    if (n1 <= n2)

        error('Index profile is anti-guiding (n1 <= n2)');

    end

    % Calls bisect - check if this is present on the path!
    if (exist('bisect','file') ~= 2)

        % Folder containing bisect.m
        addpath([userpath '/Overlap-Factor-Model/Numeric']); 

    end

    maxb = n1*k0;
    minb = n2*k0;
    npts = 10001;

    db = (maxb - minb)/(npts - 1);
    beta = minb:db:maxb;

    % Construct array of intervals for even solutions
    %
    % feven returns the function f(beta) = tan(kappa*w/2) - gamma/kappa, 
    % where
    %
    %   gamma = sqrt(beta^2- n2^2*k0^2)
    %
    %   kappa = sqrt(n1^2*k0^2 - beta^2)
    %
    % The even solutions are given for beta when f(beta) = 0.
    %
    % NOTE that f(beta) decreases monotonically between singularities, 
    % where it flips from -inf to +inf
    fev = feven(beta, k0, n1, n2, w);

    % Start with first element of fev...
    fev1 = fev(1);

    % ... see if it is less than zero
    if (fev1 < 0)

        % If it is, we are in a region where f(beta) is negative
        x1 = -1;    

    else

        % Otherwise, f(beta) is positive
        x1 = 1;     
    end

    neven = 0;                 % Number of zeros
    evenints = zeros(npts, 2); % Array for intervals

    for n = 2:npts

        % Get next element of fev
        fevn = fev(n);

        if (fevn < 0)

            x2 = -1;    % In a region where f(beta) is negative

        else

            x2 = 1;     % In a region where f(beta) is positive

        end

        % Test if we have moved from a positive region to a negative one
        % Note that f(beta) decreases monotonically except at singular 
        % points where it flips from -inf to +inf.
        if ((x1 == 1)&&(x2 == -1))

            % If so, increment the number of solutions...
            neven = neven + 1;  

            evenints(neven, 1) = beta(n-1); % ... and store interval
            evenints(neven, 2) = beta(n);

        end

        % Set first point to previous second point
        x1 = x2;

    end

    if (neven == 0)

        error('No even solutions found');

    end

    evensols = zeros(1, neven);

    % Solve for beta on each of the solution intervals
    for n = 1:neven

        % Even solution
        beta1 = evenints(n, 1);
        beta2 = evenints(n, 2);
        % Anonymous function handle
        func = @(beta) feven(beta, k0, n1, n2, w); 
        [beta_even, OK] = bisect(func, beta1, beta2, precision); 

        if (OK)

            evensols(n) = beta_even;

        else

            error(['Invalid integral for ' num2str(n) ' even solution']);

        end


    end

    % Find highest value of beta (this will be the lowest even mode)
    beta_even = max(evensols);
    neff = beta_even/k0;

    % Construct solutions
    xmax = xspan/2.0;
    xmin = -xmax;
    npts = 10001;
    dx = (xmax - xmin)/(npts - 1.0);
    x = xmin:dx:xmax;

    mode = zeros(2, length(x));
    mode(1,:) = x;

    guide = (x >= -w/2.0)&(x <= w/2.0);

    psiA = evensol(x, beta_even, k0, n1, n2, w);
    % Shifted second guide
    psiB = evensol(x+d+w, beta_even, k0, n1, n2, w);    

    % Intensity
    I = psiA.*psiA;

    % Normalise intensity
    A2 = trapezoid1D(x, I);

    % Normalisation constant
    N = 1.0/sqrt(A2); 
    psiA = N*psiA;
    psiB = N*psiB;

    mode(2,:) = psiA;

    IAA = psiA.*psiA.*guide;
    GAA = trapezoid1D(x, IAA);

    IAB = psiA.*psiB.*guide;
    GAB = trapezoid1D(x, IAB);

    eta = (c*k0/(2.0*n1^2))*(n1^2 - n2^2)*GAB;

    % Analytical expression for the coupling coefficient
%     Wr = 1.26;
%     Cn = 83.6;                  % (ns^-1)
%     eta2 = Cn*exp(-2.0*Wr*d/w);  % Coupling coefficent (absolute value)
% 
%     a = w/2.0;
% 
%     W = a*gamma(beta_even, k0, n2);
%     U = a*kappa(beta_even, k0, n1);
%     V = a*k0*sqrt(n1^2 - n2^2);
% 
%     Cn2 = (c/(k0*n1^2))*(U*W/(a*V))^2/(W+1.0);
%     eta3 = Cn2*exp(-2.0*W*d/w);

    if (report)

        lw = 2;
        gw = 0.2;

        figure;
        plot(x, gw*guide, 'k--', 'LineWidth', lw);
        title('Single guide solutions');
        ylabel('\psi');
        xlabel('x');
        grid on;
        hold on;
        plot(x, psiA, 'LineWidth', lw);

        disp(' ');
        disp(['k0 = ' num2str(k0) ' microns^(-1)']);
        disp(['n1 = ' num2str(n1,9)]);
        disp(['n2 = ' num2str(n2,9)]);
        disp(['d = ' num2str(d) ' microns']);
        disp(['w = ' num2str(w) ' microns']);
        disp(' ');
        disp(['GAA = ' num2str(GAA)]);
%         disp(['GAB = ' num2str(GAB)]);
        disp(['neff = ' num2str(neff)]);
        disp(['eta = ' num2str(eta)]);
%         disp(['eta2 = ' num2str(eta2)]);
%         disp(['eta3 = ' num2str(eta3)]);
%         disp(['Wr = ' num2str(Wr)]);
%         disp(['W = ' num2str(W)]);
%         disp(['Cn = ' num2str(Cn)]);
%         disp(['Cn2 = ' num2str(Cn2)]);

    end

    % Local functions

    function r = feven(beta, k0, n1, n2, w)
        % Argument of tan function in functions to solve
        k = kappa(beta, k0, n1);
        g = gamma(beta, k0, n2);
        
        r = tan(k*w/2) - g./k;
        
    end

    function r = kappa(beta, k0, n1)
        % Effective wavevector
        r = sqrt(n1*n1*k0*k0 - beta.*beta);
        
    end

    function r = gamma(beta, k0, n2)
        % Effective wavevector
        r = sqrt(beta.*beta - n2*n2*k0*k0);
        
    end

    % Functions to construct solutions
        
    function psi = evensol(x, beta, k0, n1, n2, w)
        
        % Effective wavevectors
        k2 = n1*n1*k0*k0 - beta*beta;
        g2 = beta*beta - n2*n2*k0*k0;
        
        k = sqrt(k2);
        g = sqrt(g2);
        
        B0 = exp(-g*w/2.0)/cos(k*w/2.0);
        
        % Masks for different regions
        m1 = (x < -w/2.0);
        m2 = (x >= -w/2.0)&(x <= w/2.0);
        m3 = (x > w/2.0);
        
        psi1 = exp(g*x).*m1;
        psi2 = B0*cos(k*x).*m2;
        psi3 = exp(-g*x).*m3;
        
        psi = psi1 + psi2 + psi3;
        
    end

end

