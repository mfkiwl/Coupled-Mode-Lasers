function [tout, Nout] = compareCoupled1D(tsim, P, d, DW, param)
% compareCoupled1D Runs the temporal evolution of coupled guide dynamics for
% comparison to PRA paper
%
% Usage:
%
%   [tout, Nout] = compareCoupled1D(tsim, P, k0, n1, n2, d, w)
%
% Arguments:
%
%   tsim        integer specifying the simulation time tspan via tspan =
%               tsim/yn, where yn is the recombination time
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
%   Nout        time evolution of variables
%
%   tout        array of time values
%

% Pump parameter (P = P/Pth)
CQ = 11.4;  % Coupled mode parameter from PRA AB11351.pdf
Q = CQ*(P - 1.0) + P;
QA = Q;
QB = Q;


[tout, Nout] = runCoupled1D(tsim, QA, QB, d, DW, param, 1);
    


end






