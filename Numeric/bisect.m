function [x0,OK] = bisect(func, x1, x2, varargin)
%BISECT Solves the equation f(x) = 0 via the bisection method
%
%Syntax:
%
%   [x0, OK] = bisect(func, x1, x2, [tol]);
%
%Arguments:
%
%   func       Function handle (passed via func = @somefunction)
%   x1         Starting guess one side of zero
%   x2         Starting guess the other side of zero
%   tol        Optional tolerance on the accuracy of the solution (default
%              value is tol = 1E-7)
%
%Return values:
%
%   x0         Solution of func(x) = 0
%   OK         Boolean indicating successful calculation 
%              
%The routine will fail under the following conditions: 
%
%   (1) func is not a function handle
%   (2) [x1, x2] is an invalid integral (i.e. func(x1) and func(x2) are not
%       either side of zero.
%
%The routine will work regardless of the order of x1 and x2 or whether
%func(x) is increasing or decreasing on the integral.

    % Test if func is a function handle
    if (~isa(func, 'function_handle'))

        x0 = 0;
        OK = false;

        return

    end

    % Check if a tolerance has been passed
    if (nargin > 3)

        tol = varargin{1};

    else

        tol = 1E-7;

    end

    % Check if x2 > x1, if not, swap them
    if (x2 < x1)

        temp = x2;
        x2 = x1;
        x1 = temp;

    end

    f1 = func(x1);
    f2 = func(x2);

    % Test if f1 and f2 are either side of zero

    if (((f1 > 0) && (f2 > 0)) || ((f1 < 0) && (f2 < 0)))

        x0 = 0;
        OK = false;

        return

    end

    OK = true;

    % Algorithm assumes increasing function, so if it is decreasing, 
    % we need to swap the signs
    if (f1 > f2)

        sn = -1.0;

    else

        sn = 1.0; 

    end

    % Initial x0
    x0 = (x1 + x2)/2.0;

    while (abs(x2 - x1) > tol)

        f0 = sn*func(x0);

        % Check if solution is exactly zero
        if (f0 == 0)

            return

        end

        % Assumes solution is negative on x1 side

        if (f0 < 0)

            x1 = x0;

        else

            x2 = x0;

        end

        x0 = (x1 + x2)/2.0;

    end

end

