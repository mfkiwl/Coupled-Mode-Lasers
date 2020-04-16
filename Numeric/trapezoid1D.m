function [I,OK] = trapezoid1D(x,y)
%TRAPEZOID1D Performs numerical integration of a function defined on a fixed and possibly non-uniform grid using the trapezoidal rule
%
%Usage:
%
%   [I,OK] = trapezoid1D(x,y);
%
%Arguments:
%
%   x   1D array for independent variable
%   y   1D array containing functional values y = y(x)
%
%Returns:
%
%   I   Integral of function
%   OK  Boolean indicating successful integration (true) or an error
%       condition (false)

    OK = false;
    I = 0;

    szx = size(x);
    szy = size(y);

    if (length(szx) > 2)

        disp('Independant array is not 1D');

        return

    end

    if ((szx(1) > 1)&&(szx(2) > 1))

        disp('Independant array is not 1D');

        return

    end


    if length(szy) > 2

        disp('Dependant array is not 1D');

        return

    end

    if ((szy(1) > 1)&&(szy(2) > 1))

        disp('Dependant array is not 1D');

        return

    end

    if (length(x) ~= length(y))

        disp('Arrays are not equal length')

    end

    x1 = x(1);
    y1 = y(1);

    for n = 2:length(x)

        x2 = x(n);
        y2 = y(n);
        dx = x2 - x1;

        dI = (y1 + y2)*dx/2.0;
        I = I + dI;

        x1 = x2;
        y1 = y2;

    end
    
end

