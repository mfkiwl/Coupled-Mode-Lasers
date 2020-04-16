function [r] = simpson1D(arrayin, dx)
%SIMPSON1D Integrates the input array via Simpson's method
%   
%Syntax:
%
%   r = simpson1D(arrayin, dx);
%
%where
%
%   arrayin is the input array to integrate
%   dx is the grid size
%
%Code relying on this routine should include the following lines to
%include the 'Numeric' folder on the path:
%
%if (exist('simpson1D','file') ~= 2)
%    
%    addpath([userpath '/Numeric']); % Folder containing simpson1D.m
%    
%end

    if ~isvector(arrayin)
        error('Input array must be 1D');
    end

    if (dx < 0)
        error('dx must be greater than zero');
    end

    n = length(arrayin);

    index = 1:n;
    % Need to make sure index and arrayin have the same shape
    index = reshape(index, size(arrayin)); 
    modd = 2.0*(mod(index,2) == 1);
    meven = 4.0*(mod(index,2) == 0);
    modd(1) = 1.0;
    modd(end) = 1.0;
    mask = modd + meven;

    r = sum(mask.*arrayin)*dx/3.0;

end

