% ورودی‌ها
a = 40;
b = 40;
c = 30;
R1 = 30;
R2 = 30;
t = 10;
sigma = 324;  % MPa
ro = 2780;    % kg/m³
E = 73;    % GPa

% محاسبات اولیه
H = 2*t + R1 + R2 + a;
Yc = ((4*R2 + 4*R1 + 2*b + 2*a + 8*t)*(a + R1 + R2 + 2*t) - c*(1.5*t + a + R1 + R2)) / ...
     (8*t + 4*(R1 + R2) + 2*(a + b) - c);
I = (2*b*t^3 + 2*t*a^3 - c*t^3)/12 - c*t*(Yc - (R1 + R2 + a + 1.5*t))^2 + ...
    b*t*((Yc - 0.5*t)^2 + (Yc - (a + R1 + R2 + 1.5*t))^2) + ...
    2*a*t*(Yc - (R1 + R2 + 2*t + a)/2)^2;

if (H - Yc) < Yc
    C = Yc;
else
    C = H - Yc;
end

A = t*(2*b + 2*a - c + 3.14159*(R1 + R2 + t));
F = C/I;
Q = 4*F*10^3*(500 + 1.8*10^(-6)*ro*A);
Y = 10^6*(2000/48 + 100*ro*A*(0.0002 + 5/384))/(E*I);

fprintf('F = %.4f\n', F);
fprintf('Q = %.4f MPa\n', Q);
fprintf('Y = %.4f mm\n', Y);

% بهینه‌سازی
x0 = [a, b, c, R1, R2, t];
lb = [5, 5, 30, 3, 3, 1];
ub = [10000, 10000, 10000, 1000, 1000, 10];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 500);

[x_opt, fval, exitflag] = fmincon(@(x) objective_func(x, sigma, ro, E), x0, ...
    [], [], [], [], lb, ub, @(x) nonlinear_constraints(x, sigma, ro, E), options);

if exitflag > 0
    [A_opt, Q_opt, Y_opt] = evaluate_func(x_opt, sigma, ro, E);
    fprintf('\n--- طرح بهینه (حداقل مساحت) ---\n');
    fprintf('a = %.2f, b = %.2f, c = %.2f\n', x_opt(1), x_opt(2), x_opt(3));
    fprintf('R1 = %.2f, R2 = %.2f, t = %.2f\n', x_opt(4), x_opt(5), x_opt(6));
    fprintf('مساحت = %.2f mm²\n', A_opt);
    fprintf('Q = %.2f MPa  (باید ≤ %.2f باشد)\n', Q_opt, sigma);
    fprintf('Y = %.2f mm  (باید ≤ 5 باشد)\n', Y_opt);
else
    fprintf('بهینه‌سازی موفق نبود\n');
end

% توابع
function [A, Q, Y] = evaluate_func(vars, sigma, ro, E)
    a = vars(1); b = vars(2); c = vars(3);
    R1 = vars(4); R2 = vars(5); t = vars(6);
    
    H = 2*t + R1 + R2 + a;
    Yc = ((4*R2 + 4*R1 + 2*b + 2*a + 8*t)*(a + R1 + R2 + 2*t) - c*(1.5*t + a + R1 + R2)) / ...
         (8*t + 4*(R1 + R2) + 2*(a + b) - c);
    I = (2*b*t^3 + 2*t*a^3 - c*t^3)/12 - c*t*(Yc - (R1 + R2 + a + 1.5*t))^2 + ...
        b*t*((Yc - 0.5*t)^2 + (Yc - (a + R1 + R2 + 1.5*t))^2) + ...
        2*a*t*(Yc - (R1 + R2 + 2*t + a)/2)^2;
    
    if (H - Yc) < Yc
        C = Yc;
    else
        C = H - Yc;
    end
    
    A = t*(2*b + 2*a - c + 3.14159*(R1 + R2 + t));
    F = C/I;
    Q = 4*F*10^3*(500 + 1.8*10^(-6)*ro*A);
    Y = 10^6*(2000/48 + 100*ro*A*(0.0002 + 5/384))/(E*I);
end

function obj = objective_func(vars, sigma, ro, E)
    [A, ~, ~] = evaluate_func(vars, sigma, ro, E);
    obj = A;
end

function [c_ineq, c_eq] = nonlinear_constraints(vars, sigma, ro, E)
    [~, Q, Y] = evaluate_func(vars, sigma, ro, E);
    a = vars(1); b = vars(2); c_var = vars(3);
    R1 = vars(4); R2 = vars(5); t = vars(6);
    
    % قیود نامساوی (باید ≤ 0)
    c_ineq = [
        Q - sigma;           % تنش
        Y - 5;               % خیز
        30 - c_var;          % c حداقل
        3*t - R1;            % R1 حداقل
        3*t - R2;            % R2 حداقل
        c_var - (b + 2*R1);  % c حداکثر
    ];
    
    c_eq = [];
end