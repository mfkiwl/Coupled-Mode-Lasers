function integral = trapezoid2D(x, y, fxy)
%TRAPEZOID2D Uses a generalization of the trapezoidal method to integrate a 2D grid with a possible nonequal grid size 
%  
%Usage:
%
%   integral = trapezoid2D(x, y, fxy);
%
%Arguments:
%
%   x       independent x variable
%
%   y       independent y variable
%
%   fxy     2D function f(x,y)
%
%Notes:
%
%The input array is assumed to be in matrix order, so that the first index
%(row index) corresponds to the y coordinate and the second index (column
%index) to the x coordinate.
%
%That is:
%
%   fxy(1,1)     = f(xmin, ymin)
%   fxy(end,1)   = f(xmin, ymax)
%   fxy(1,end)   = f(xmax, ymin)
%   fxy(end,end) = f(xmax, ymax)

    % Check dimensions
    xlen = length(x);
    ylen = length(y);
    sz = size(fxy);
    m = sz(1);
    n = sz(2);

    if ((xlen ~= n)||(ylen ~= m))

        error('Function array dimensions do not match x and y arrays');

    end

    integral = 0.0;


    for j = 1:n-1

        dx = x(j+1) - x(j);

        for i = 1:m-1

            dy = y(i+1) - y(i);

            dA = dx*dy;

            fdA = (fxy(i+1,j+1) + fxy(i+1,j) + fxy(i,j+1) + fxy(i,j))*dA/4.0;

            integral = integral + fdA;

        end

    end
        
end

