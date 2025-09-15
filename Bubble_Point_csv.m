%% Bubble-point calculation using CSV for components
clc; clear;

%% --- Load component data ---
components_table = readtable('data/component.csv');  % read CSV

% --- Select which components to use (choose row numbers) ---
selected = [1 2];   % example: n-butane + n-pentane

names = components_table.Name(selected);
Tc = components_table.Tc(selected);
Pc = components_table.Pc(selected);
omega = components_table.omega(selected);

n = length(selected);  % number of components

%% --- User-defined conditions ---
T = 350;              % temperature [K]
x = [0.5, 0.5];       % liquid mole fractions
P_guess = 5e5;        % initial pressure [Pa]

%% --- Iterative bubble-point calculation ---
P = P_guess;
tol = 1e-5;
max_iter = 100;

for iter = 1:max_iter
    % Step 1: calculate a and b for each component
    a = zeros(1,n);
    b = zeros(1,n);
    for i = 1:n
        [a(i), b(i)] = pr_params(T, Tc(i), Pc(i), omega(i));
    end
    
    % Step 2: calculate Z for vapor and liquid
    Zv = zeros(1,n); Zl = zeros(1,n);
    for i = 1:n
        Zroots = pr_Z_roots(P, T, a(i), b(i));
        Zv(i) = max(Zroots);   % vapor root
        Zl(i) = min(Zroots);   % liquid root
    end
    
    % Step 3: calculate fugacity coefficients
    phi_v = zeros(1,n); phi_l = zeros(1,n);
    for i = 1:n
        phi_v(i) = exp(pr_lnphi(P, T, a(i), b(i), Zv(i)));
        phi_l(i) = exp(pr_lnphi(P, T, a(i), b(i), Zl(i)));
    end
    
    % Step 4: calculate y_i
    y = x .* phi_l ./ phi_v;
    y = y / sum(y);  % normalize

    % Step 5: update pressure
    P_new = P * sum(x .* phi_l ./ phi_v);
    
    % Convergence check
    if abs(P_new - P)/P < tol
        P = P_new;
        break
    else
        P = P_new;
    end
end

%% --- Display results ---
disp(['Bubble-point pressure = ', num2str(P/1e5), ' bar']);
disp('Vapor composition y_i = ');
for i = 1:n
    disp([names{i}, ': ', num2str(y(i))]);
end
