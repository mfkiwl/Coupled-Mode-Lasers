function [neff, modes] = singleSlab(guide)
%SINGLESLAB Calculates the effective indices and modal profiles for a
%single real slab waveguide.
% 
%%  Description
%
%
%%  Usage
%
%    [neff, modes] = singleSlab(guide);
%
%%  Arguments
%
%    guide      A structure containing the waveguide parameters:
%
%               guide.k0        free space wavevector (1/micron)
%               guide.ncore     core refractive index
%               guide.nclad     cladding refractive index
%               guide.width     guide width (micron)
%               guide.span      total spatial span to calculate modal
%                               solutions on (typically 7*width)
%
%%  Returns
%
%   neff        (ntot x 1) array of the effective indices of the solutions,
%               where ntot is the total number of confined solutions found
%               (lowest modal solution first)
%
%   modes       (ntot x npts) spatial array containing the modal profiles 
%               (npts is currently hard-coded to 10001)
%
%%  Notes
%
%
%%  Code

    k0 = guide.k0;
    n1 = guide.ncore;
    n2 = guide.nclad;
    w = guide.width;
    xspan = guide.span;

    precision = 1E-9;

    % If profile is anti-guiding, report an error
    if (n1 <= n2)

        error('Index profile is anti-guiding (n1 <= n2)');

    end
    

    % Calls bisect - check if this is present on the path!
    if (exist('bisect','file') ~= 2)

        addpath([userpath '/Numeric']); % Folder containing bisect.m

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

        x1 = -1;    % If it is, we are in a region where f(beta) is negative

    else

        x1 = 1;     % Otherwise, f(beta) is positive
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

            neven = neven + 1;  % If so, increment the number of solutions...

            evenints(neven, 1) = beta(n-1); % ... and store interval
            evenints(neven, 2) = beta(n);

        end

        % Set first point to previous second point
        x1 = x2;

    end

    if (report)
        if (neven == 0)

            disp('No even solutions found');

        else

            if (neven ==1)
                estr = ' even solution';
            else
                estr = ' even solutions';
            end

            disp(['Found ' num2str(neven) estr]);

        end
    end

    % Construct array of intervals for odd solutions
    %
    % fod returns the function f(beta) = cot(kappa*w/2) + gamma/kappa, 
    % where
    %
    %   gamma = sqrt(beta^2- n2^2*k0^2)
    %
    %   kappa = sqrt(n1^2*k0^2 - beta^2)
    %
    % The odd solutions are given for beta when f(beta) = 0.
    %
    % NOTE that f(beta) increases monotonically between singularities, 
    % where it flips from +inf to -inf
    fod = fodd(beta, k0, n1, n2, w);

    if (fod(1) < 0)
        x1 = -1;
    else
        x1 = 1;
    end

    nodd = 0;                 % Number of zeros
    oddints = zeros(npts, 2); % Array for intervals

    for n = 2:npts

        if (fod(n) < 0)
            x2 = -1;
        else
            x2 = 1;
        end

        if ((x1 == -1)&&(x2 == 1))

            nodd = nodd + 1;
            % Store interval
            oddints(nodd, 1) = beta(n-1);
            oddints(nodd, 2) = beta(n);   

        end

        x1 = x2;

    end

    evensols = zeros(1, neven);
    oddsols = zeros(1, nodd);

    % Solve for beta on each of the solution intervals
    for n = 1:neven

        % Even solution
        beta1 = evenints(n, 1);
        beta2 = evenints(n, 2);
        func = @(beta) feven(beta, k0, n1, n2, w); % Anonymous function handle
        [beta_even, OK] = bisect(func, beta1, beta2, precision); 

        if (OK)

            evensols(n) = beta_even;

        else

            warning(['Invalid integral for ' num2str(n) ' even solution']);

        end


    end

    for n = 1:nodd

        % Odd solution
        beta1 = oddints(n, 1);
        beta2 = oddints(n, 2);
        func = @(beta) fodd(beta, k0, n1, n2, w); % Anonymous function handle
        [beta_odd, OK] = bisect(func, beta1, beta2, precision); 

        if (OK)

            oddsols(n) = beta_odd;

        else

            warning(['Invalid integral for ' num2str(n) ' odd solution']);

        end


    end

    % Compile total solutions
    ntot = neven + nodd;
    neff = zeros(ntot, 1);

    for n = 1:ntot

        if (mod(n,2) == 1)

            % Store even solutions
            m = (n+1)/2;
            neff(n) = evensols(m)/k0;

        else

            % Store odd solutions
            m = n/2;
            neff(n) = oddsols(m)/k0;

        end

    end

    % Construct solutions
    xmax = xspan/2.0;
    xmin = -xmax;
    npts = 10001;
    dx = (xmax - xmin)/(npts - 1.0);
    x = xmin:dx:xmax;

    modes = zeros(ntot+1, length(x));
    modes(1,:) = x;

    % guide = (x >= -w/2.0)&(x <= w/2.0);

    for m = 1:neven

        n = 2*m - 1;

        beta_even = evensols(m);
        psi_s = evensol(x, beta_even, k0, n1, n2, w);

        % Intensity
        I = psi_s.*psi_s;

        % Normalise intensity
        A2 = trapezoid1D(x, I);

        % Normalisation constant
        N = 1.0/sqrt(A2); 
        psi_s = N*psi_s;

        modes(n+1,:) = psi_s;

    end

    for m = 1:nodd

        n = 2*m;

        beta_odd = oddsols(m);
        psi_a = oddsol(x, beta_odd, k0, n1, n2, w);

        % Intensity
        I = psi_a.*psi_a;

        % Normalise intensity
        A2 = trapezoid1D(x, I);

        % Normalisation constant
        N = 1.0/sqrt(A2); 
        psi_a = N*psi_a;

        modes(n+1,:) = psi_a;

    end
    
    % --------------------- FUNCTION RETURNS ------------------------------

    
    % ---------------------------------------------------------------------
    %       LOCAL FUNCTIONS
    % ---------------------------------------------------------------------

    % Used to find even solutions
    function r = feven(beta, k0, n1, n2, w)
        % Argument of tan function in functions to solve
        k = kappa(beta, k0, n1);
        g = gamma(beta, k0, n2);
        
        r = tan(k*w/2) - g./k;
        
    end

    % Used to find odd solutions
    function r = fodd(beta, k0, n1, n2, w)
        % Argument of tan function in functions to solve
        k = kappa(beta, k0, n1);
        g = gamma(beta, k0, n2);
        
        r = cot(k*w/2) + g./k;
        
    end

    function r = kappa(beta, k0, n1)
        % Effective wavevector
        r = sqrt(n1*n1*k0*k0 - beta.*beta);
        
    end

    function r = gamma(beta, k0, n2)
        % Effective wavevector
        r = sqrt(beta.*beta - n2*n2*k0*k0);
        
    end

    % ---------------------------------------------------------------------
    %            FUNCTIONS TO CONSTRUCT MODAL PROFILES
    % ---------------------------------------------------------------------
    
    % Even solution
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

    % Odd solution
    function psi = oddsol(x, beta, k0, n1, n2, w)
        
        % Effective wavevectors
        k2 = n1*n1*k0*k0 - beta*beta;
        g2 = beta*beta - n2*n2*k0*k0;
        
        k = sqrt(k2);
        g = sqrt(g2);
        
        B0 = -exp(-g*w/2.0)/sin(k*w/2.0);
        
        % Masks for different regions
        m1 = (x < -w/2.0);
        m2 = (x >= -w/2.0)&(x <= w/2.0);
        m3 = (x > w/2.0);
        
        psi1 = exp(g*x).*m1;
        psi2 = B0*sin(k*x).*m2;
        psi3 = -exp(-g*x).*m3;
        
        psi = psi1 + psi2 + psi3;
        
    end

end

