function [r1, r2] = gaussian(nsize)
%GAUSSIAN Summary of this function goes here
%   

    r1 = randn(nsize,1);
    r2 = randn(nsize,1);
    
    figure;
    h = histogram(r1);
    hold on;

    npts = 201;
    xmin = h.BinLimits(1);
    xmax = h.BinLimits(2);
    dx = (xmax - xmin)/(npts-1);
    x = xmin:dx:xmax;

    dr = h.BinWidth;

    y = nsize*dr*exp(-x.*x/2)/sqrt(2*pi);

    plot(x, y, 'LineWidth', 2);
    
    figure;
    histogram(r2);
    hold on;
    plot(x, y, 'LineWidth', 2);
    
    var = sum(r1.*r1)/nsize;
    covar = sum(r1.*r2)/nsize;
    
    disp(['var = ' num2str(var)])
    disp(['covar = ' num2str(covar)])
    





end

