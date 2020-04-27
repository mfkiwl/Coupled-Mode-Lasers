function [J] = besselj_test(nu)
%BESSELJ_TEST Testing routine for Bessel functions
%   

    % Domain for Bessel function
    z0 = 0.0;
    z1 = 20;
    npts = 201;
    dz = (z1 - z0)/(npts - 1);
    
    z = z0:dz:z1;

    J = besselj(nu,z);
    
    lw = 1.5;   % Linewidth
    
    figure
    plot(z,J,'LineWidth', lw)
    grid on
    
    title(['Bessel Function of the First Kind for $\nu = ' num2str(nu) '$'],'interpreter','latex')
    xlabel('$z$','interpreter','latex')
    ylabel('$J_\nu(z)$','interpreter','latex')
    
    % Derivative of the Bessel function
    dJ_num = (J(3:npts) - J(1:(npts-2)))/(2*dz);
    z2 = z(2:(npts-1));
    
    dJ = deriv_besselj(nu, z);
    
    figure
    grid on
    hold on
    plot(z2,dJ_num, 'd')
    plot(z,dJ,'LineWidth', lw)

    title(['Derivative of Bessel Function of the First Kind for $\nu = ' num2str(nu) '$'],'interpreter','latex')
    xlabel('$z$','interpreter','latex')
    ylabel('$dJ_\nu(z)/dz$','interpreter','latex')
    
    legend('Numerical', 'Analytic')
    
    


end

