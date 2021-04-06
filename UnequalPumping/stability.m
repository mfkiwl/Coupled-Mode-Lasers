
% Calls bisect - check if this is present on the path!
if (exist('bisect','file') ~= 2)

    % Folder containing bisect.m
    addpath([userpath '/Coupled Mode Lasers/Numeric']); 

end


QBmax = 25;
QBmin = 1.125;
npts = 201;
dQ = (QBmax - QBmin)/(npts - 1);

QB = QBmin:dQ:QBmax;

QA = 3.2;

q = (QA - QB)./(QA + QB - 2);

phis = zeros(size(q));

eta = 1;
alpha = 4;
tau_N = 2;
tau_p = 0.002;

lw = 2;
figure
hold on

phi_max = pi;
phi_min = -pi;
ppts = 21;

dphi = (phi_max - phi_min)/(ppts - 1);

phi = phi_min:dphi:phi_max;

for n=1:ppts
    
    plot(q, eq(phi(n), q, QA, QB, eta, alpha, tau_N), 'LineWidth', lw)
    
    drawnow
    
    pause
    
    
%     % Anonymous function handle
%     func = @(phi) eq(phi, q(n), QA, QB(n), eta, alpha, tau_N); 
%     
%     [ph, OK] = bisect(func, -3*pi/2, 3*pi/2);
%     
%     if (OK)
%         phis(n) = ph;
%     else
%         phis(n) = NaN;
%     end

    
end

grid on


% DW = 2*eta*sqrt((alpha*alpha + q.*q)./(1 - q.*q));
% 
% 
% figure
% hold on
% plot(real(DW), QB, 'k', 'LineWidth', lw);
% plot(-real(DW), QB, 'k', 'LineWidth', lw);
% xlim([-15 15])


function r = eq(phi, q, QA, QB, eta, alpha, tau_N)

    r = ((QA+QB).*(1-q.*q) + 2*q.*q).*sqrt(1-q.*q) + 4*eta*tau_N*(alpha*cos(phi) + q*sin(phi));
end