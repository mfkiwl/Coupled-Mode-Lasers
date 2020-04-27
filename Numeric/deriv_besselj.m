function [dJ] = deriv_besselj(nu, z)
%DERIV_BESSELJ Returns the derivative of the Bessel function J(nu,z)
%%
% *DERIV_BESSELJ*
%
%%  Description
%
% Returns the derivative of the Bessel function J(nu, z)
%
% From Abramowitz and Stegun, p. 361, 9.1.27, the Bessel functions 
% J(α), Y(α), H(α)(1) and H(α)(2) all satisfy the recurrence relations
%
%   dZ(α,x)/dx = (1/2)(Z(α-1,x) - Z(α+1,x))
%
%%  Usage
%
%   [dJ] = deriv_besselj(nu, z)
%
%%  Arguments
%
%   nu      the order of the Bessel function
%   z       array containing the domain of the function
%
%%  Returns
%
%   dJ      array containing dJ(nu, z)/dz on the domain z
%
%%  Code

    dJ = 0.5*(besselj(nu-1, z) - besselj(nu+1, z));
    
end

