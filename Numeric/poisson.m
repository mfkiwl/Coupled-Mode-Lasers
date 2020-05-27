function poisson(lambda)
%POISSON Summary of this function goes here
%   

    N = length(lambda);
    kmax = 100;
    kmin = 0;

    k = kmin:kmax;

    figure
    hold on

    for n = 1:N

        % Generate Poisson distribution
        lam = lambda(n);
        Pk = (lam.^k)*exp(-lam)./factorial(k);
        plot(k, Pk, '-d')

    end

