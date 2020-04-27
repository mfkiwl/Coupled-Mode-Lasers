function [dK] = deriv_besselk(nu, z)
%DERIV_BESSELK Returns the derivative of the Bessel function K(nu,z)
%%
% *DERIV_BESSELK*
%
%%  Description
%
% Returns the derivative of the Bessel function K(nu, z)
%
%
%%  Usage
%
%   [dK] = deriv_besselk(nu, z)
%
%%  Arguments
%
%   nu      the order of the Bessel function
%   z       array containing the domain of the function
%
%%  Returns
%
%   dK      array containing dK(nu, z)/dz on the domain z
%
%%  Code

    dK = -0.5*(besselk(nu-1, z) + besselk(nu+1, z));
    
end

