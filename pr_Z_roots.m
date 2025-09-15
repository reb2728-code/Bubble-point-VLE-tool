function Z = pr_Z_roots(P, T, a, b)
% pr_Z_roots calculates the compressibility factor(s) Z for Peng-Robinson EOS
% Inputs:
%   P = pressure [Pa]
%   T = temperature [K]
%   a = Peng-Robinson attraction parameter [Pa*m^6/mol^2]
%   b = Peng-Robinson co-volume parameter [m^3/mol]
% Output:
%   Z = compressibility factor (liquid and vapor)

R = 8.314462618; % Gas constant [J/(mol·K)]

% Step 1: calculate dimensionless parameters
A = a * P / (R^2 * T^2);
B = b * P / (R * T);

% Step 2: define the cubic equation coefficients
% Z^3 - (1-B)*Z^2 + (A - 3*B^2 - 2*B)*Z - (A*B - B^2 - B^3) = 0
coeffs = [1, -(1-B), A - 3*B^2 - 2*B, -(A*B - B^2 - B^3)];

% Step 3: solve the cubic equation
cubic_roots = roots(coeffs);  % roots will be found

% Step 4: keep only real roots (physical solutions)
Z = cubic_roots(imag(cubic_roots)==0); 
end
