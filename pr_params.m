function [a,b] = pr_params(T, Tc, Pc, omega)
% pr_params calculates Peng-Robinson EOS parameters for a pure component
% Inputs:
%   T    = temperature [K]
%   Tc   = critical temperature [K]
%   Pc   = critical pressure [Pa]
%   omega = acentric factor
% Outputs:
%   a = attraction parameter [Pa*m^6/mol^2]
%   b = co-volume parameter [m^3/mol]

R = 8.314462618; % J/(mol·K)

% Step 1: calculate kappa
kappa = 0.37464 + 1.54226*omega - 0.26992*omega^2;

% Step 2: calculate alpha(T)
alpha = (1 + kappa*(1 - sqrt(T/Tc)))^2;

% Step 3: calculate 'a' and 'b'
a = 0.45724 * (R^2 * Tc^2 / Pc) * alpha;
b = 0.07780 * (R * Tc / Pc);
end
