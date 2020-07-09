function [eta] = realEta(d_over_a, a, k0, n1, n2)
%REALETA Returns the real part of the coupling coefficient in the coupled
%mode model as function of d/a, the ratio of guide separation to guide
%width for slab guides
%  
%% Usage
%
%   [eta] = realEta(d_over_a, a, k0, n1, n2);
%
%% Arguments
%
%   d_over_a    an array containing the values of d/a to calculate for
%   a           guide width
%   k0          free space wavevector (1/micron)
%   n1          core refractive index
%   n2          cladding refractive index
%
%% Returns
%
%   eta         an array (size of da) containing the real coupling
%               coefficient
%
%% Notes:
%
%   Calls singleSlab(k0, n1, n2, w, d)
%

    eta = zeros(size(d_over_a));
    
    for n = 1:length(eta)
        
        d = a*d_over_a(n);
        [~, eta_n, ~, ~] = singleSlab(k0, n1, n2, a, d);
        
        eta(n) = eta_n;
        
    end
     
end

