function [beta_even, Gamma] = slabGuideBeta(k0, n1, n2, w)
%SLABGUIDEBETA Calculates parameters for a slab waveguide.
%%
% *SLABGUIDEBETA*
%
%%  Description
%
% This routine calculates the propagation constant for a slab guide.
% Only the lowest even solution found is returned. 
% 
%%  Usage
%
%   [beta] = slabGuideBeta(k0, n1, n2, w);
%
%%  Arguments
%
%   k0          free space wavevector (1/micron)
%   n1          core refractive index
%   n2          cladding refractive index
%   w           guide width (micron)
%
%%  Returns
%
%   beta_even   the propagation constant for the guide for wavevector k0
%   Gamma       optical confinement factor
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
%   to add it. 
%
%% Throws
%
%   Possible exceptions thrown by this routine: 
%
%       Index profile is anti-guiding (ncore <= nclad)
%
%       No even solutions found
%
%       Invalid integral for <n> even solution (error in bisection routine)
%
%% Code

    % precision used in bisection routine
    precision = 1E-9;

    % If profile is anti-guiding, report an error
    if (n1 <= n2)

        errID = 'Bisect:Exception';
        errMsg = 'Index profile is anti-guiding (ncore <= nclad)';
        ME = MException(errID, errMsg);
        throw(ME);

    end

    % Calls bisect - check if this is present on the path!
    if (exist('bisect','file') ~= 2)

        % Folder containing bisect.m
        addpath([userpath '/Overlap-Factor-Model/Numeric']); 

    end

    % Set up array of values for propagation constant beta
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

        errID = 'Bisect:Exception';
        errMsg = 'No even solutions found';
        ME = MException(errID, errMsg);
        throw(ME);

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
            
            errID = 'Bisect:Exception';
            errMsg = 'Invalid interval';
            ME = MException(errID, errMsg);
            throw(ME);
        
        end


    end

    % Find highest value of beta (this will be the lowest even mode)
    beta_even = max(evensols);
    
    % Construct solutions
    xspan = 9.0*w;
    xmax = xspan/2.0;
    xmin = -xmax;
    npts = 10001;
    dx = (xmax - xmin)/(npts - 1.0);
    x = xmin:dx:xmax;

    guide = (x >= -w/2.0)&(x <= w/2.0);

    psiA = evensol(x, beta_even, k0, n1, n2, w);
    
    % Intensity
    I = psiA.*psiA;

    % Normalise intensity
    A2 = trapezoid1D(x, I);

    % Normalisation constant
    N = 1.0/sqrt(A2); 
    psiA = N*psiA;

    IAA = psiA.*psiA.*guide;
    Gamma = trapezoid1D(x, IAA);
    
    
%% Local functions

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

