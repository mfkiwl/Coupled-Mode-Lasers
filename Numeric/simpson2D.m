function [r] = simpson2D(arrayin, dx, dy)
%SIMPSON2D Integrates a 2D input array via Simpson's method
%   
%Syntax:
%
%   r = simpson2D(arrayin, dx, dy);
%
%where
%
%   arrayin is the input array to integrate 
%   dx is the grid size in the x direction
%   dy is the grid size in the y direction
%
%   A row of arrayin is taken to be over the x axis (so each row labels a
%   different y coordinate). That is, the array has the form z(row, col)
%   where
%
%   z(x1, y1) z(x2, y1) ... z(xm, y1)
%   z(x1, y2) z(x2, y2) ... z(xm, y2)
%       .
%       .
%       .
%   z(x1, yn) z(x2, yn) ... z(xm, yn)
%
%Note: Calls simpson1D to perform each 1D integral
%
%Code relying on this routine should include the following lines to
%include the 'Numeric' folder (which includes both simpson1D and simpson2D)
%on the path:
%
%if (exist('simpson1D','file') ~= 2)
%    
%    addpath([userpath '/Numeric']); % Folder containing simpson1D.m
%    
%end
%
%Example:
%
%   dx = pi/100;
%   dy = pi/200;
%   x = 0:dx:100*dx;        % x = 0 to pi over 101 grid points
%   y = 0:dy:100*dy;        % y = 0 to pi/2 over 101 grid points
%   [X, Y] = meshgrid(x,y); % Set up 2D grid
%   z = sin(X).*sin(Y);     % Set up test function 
%   surf(z);                % Plot function to see what it looks like
%   xlabel('x');
%   ylabel('y');
%   r = simpson2D(z, dx, dy);   % Should return (and does) r = 2
%


    if (dx < 0)
        error('dx must be greater than zero');
    end

    if (dy < 0)
        error('dy must be greater than zero');
    end

    % Squeeze out any singleton dimensions
    arrayin = squeeze(arrayin);

    sz = size(arrayin);

    if ((length(sz) > 2)||isvector(arrayin))

        error('Input array must be 2D');

    end

    rows = sz(1);

    % Create a new 1D array
    array1D = zeros(rows, 1);

    for n = 1:rows

        % Integrate each row over x
        array1D(n) = simpson1D(arrayin(n,:), dx);

    end

    % Integrate the rows over y
    r = simpson1D(array1D, dy);

end

