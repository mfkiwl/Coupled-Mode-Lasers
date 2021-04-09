function [param] = loadParams(filename)
%LOADPARAMS Loads parameters from MAT file
%%
% *LOADPARAMS*
%
%%  Description
%
% Loads parameters from a MAT file and returns them in a structure for use
% with other routines.
%  
%%  Usage
%
%   [param] = loadParams(filename)
%
%%  Arguments
%
%   filename    the filename of the MAT file as a string
%
%%  Returns
%
%   param	a structure containing the loaded parameters. The initialised
%           fields are:
%               param.aH	linewidth enhancement factor 
%               param.k0    free-space wavevector (1/micron)
%               param.kp    cavity loss rate 1/(2*tau_p) (1/ns) 
%               param.n1    refractive index in the core 
%               param.n2    refractive index in the cladding 
%               param.w     width of waveguide (microns) 
%               param.yn    inverse carrier lifetime 1/tau_N (1/ns)
%%

    param = load(filename);


end

