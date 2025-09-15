% Components: n-butane and n-pentane
components = {'n-butane','n-pentane'};
Tc = [425.2, 469.7];   % K
Pc = [3796000, 3370000]; % Pa
omega = [0.2, 0.251];

% Liquid mole fractions (x_i)
x = [0.5, 0.5]; % assume equimolar mixture

% Temperature [K] for bubble-point
T = 350;

% Initial guess for pressure [Pa]
P_guess = 5e5;

P = P_guess;        % initial pressure
tol = 1e-5;         % convergence tolerance
max_iter = 100;     % maximum iterations

for iter = 1:max_iter
    % Step 1: calculate a and b for each component
    for i = 1:2
        [a(i), b(i)] = pr_params(T, Tc(i), Pc(i), omega(i));
    end
    
    % Step 2: calculate Z for vapor (largest root) and liquid (smallest root)
    Zv = zeros(1,2);  % vapor Z for each component
    Zl = zeros(1,2);  % liquid Z
    for i = 1:2
        Zroots = pr_Z_roots(P, T, a(i), b(i));
        Zv(i) = max(Zroots); % vapor root
        Zl(i) = min(Zroots); % liquid root
    end
    
    % Step 3: calculate fugacity coefficients
    phi_v = zeros(1,2);
    phi_l = zeros(1,2);
    for i = 1:2
        phi_v(i) = exp(pr_lnphi(P, T, a(i), b(i), Zv(i))); % convert ln(phi) to phi
        phi_l(i) = exp(pr_lnphi(P, T, a(i), b(i), Zl(i)));
    end
    
    % Step 4: calculate y_i using phi_i
    y = x .* phi_l ./ phi_v;
    y = y / sum(y); % normalize mole fractions
    
    % Step 5: update pressure estimate
    P_new = P * sum(x .* phi_l ./ phi_v); 
    
    % Check convergence
    if abs(P_new - P)/P < tol
        P = P_new;
        break
    else
        P = P_new;
    end
end

disp(['Bubble-point pressure = ', num2str(P/1e5), ' bar']);
disp('Vapor composition y_i = ');
disp(y);
