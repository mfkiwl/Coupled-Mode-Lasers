function [J] = bessel_test(nu)
%BESSEL_TEST Testing routine for Bessel functions
%   

    % Domain for Bessel function
    z0 = 0.0;
    z1 = 20;
    npts = 201;
    dz = (z1 - z0)/(npts - 1);
    
    z = z0:dz:z1;

    J = zeros(1, npts);
    
    J(1,:) = besselj(nu,z);
    
    figure
    plot(z,J)
    grid on

    title(['Bessel Function of the First Kind for $\nu = ' num2str(nu) '$'],'interpreter','latex')
    xlabel('$z$','interpreter','latex')
    ylabel('$J_\nu(z)$','interpreter','latex')
    
    


end

