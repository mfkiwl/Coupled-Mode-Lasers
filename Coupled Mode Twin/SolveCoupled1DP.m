function [Ns, found, E, esign] = SolveCoupled1DP(P, d, DW, param, varargin)
%SolveCoupled1DP Finds steady state solutions for coupled mode model
%
% Usage:
%
%   [Ns, E, found] = SolveCoupled1DP2(P, d, DW, param, varargin)
%
% Arguments:
%
%   P           P/Pth pump parameter from coupled mode model
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

% Pump parameter (P = P/Pth)
CQ = 11.4;  % Coupled mode parameter from PRA AB11351.pdf
Q = CQ*(P - 1.0) + P;
QA = Q;
QB = Q;

if (nargin >  4)
    
    [Ns, found, E, esign] = SolveCoupled1D(QA, QB, d, DW, param, 1);
    
    
else
    

    [Ns, found, E, esign] = SolveCoupled1D(QA, QB, d, DW, param);
    
end
    


end






