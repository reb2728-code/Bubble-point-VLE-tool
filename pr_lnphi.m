function lnphi = pr_lnphi(P, T, a, b, Z)
% pr_lnphi calculates ln(fugacity coefficient) for a pure component
% Inputs:
%   P = pressure [Pa]
%   T = temperature [K]
%   a = PR attraction parameter [Pa*m^6/mol^2]
%   b = PR co-volume parameter [m^3/mol]
%   Z = compressibility factor (from pr_Z_roots)
% Output:
%   lnphi = natural log of fugacity coefficient

R = 8.314462618; % J/(mol·K)

A = a * P / (R^2 * T^2);
B = b * P / (R * T);

% PR EOS fugacity equation for a pure component
lnphi = Z - 1 - log(Z - B) - (A/(2*sqrt(2)*B))*log((Z + (1+sqrt(2))*B)/(Z + (1-sqrt(2))*B));

end
