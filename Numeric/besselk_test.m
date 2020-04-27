function [K] = besselk_test(nu)
%BESSELK_TEST Testing routine for Bessel functions
%   

    % Domain for Bessel function
    z0 = 0.0;
    z1 = 1.0;
    npts = 201;
    dz = (z1 - z0)/(npts - 1);
    
    z = z0:dz:z1;

    K = besselk(nu,z);
    
    lw = 1.5;   % Linewidth
    
    figure
    plot(z,K,'LineWidth', lw)
    grid on
    
    title(['Bessel Function $K_{' num2str(nu) '}$'],'interpreter','latex')
    xlabel('$z$','interpreter','latex')
    ylabel('$K_\nu(z)$','interpreter','latex')
    
    %ylim([0 20])
    
    % Derivative of the Bessel function
    dK_num = (K(3:npts) - K(1:(npts-2)))/(2*dz);
    z2 = z(2:(npts-1));
    
    dK = deriv_besselk(nu, z);
     
    figure
    grid on
    hold on
    plot(z2,dK_num, 'd')
    plot(z,dK,'LineWidth', lw)

    title(['Derivative of Bessel Function $K_{' num2str(nu) '}$'],'interpreter','latex')
    xlabel('$z$','interpreter','latex')
    ylabel('$dK_\nu(z)/dz$','interpreter','latex')
    
    legend('Numerical', 'Analytic')
    
    


end

