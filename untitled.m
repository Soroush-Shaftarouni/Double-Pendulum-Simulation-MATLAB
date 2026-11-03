% ورودی‌های اولیه
material = input('1 for Al, 2 for Steel: ');
if material == 1
    sigma = 324;
    ro = 2780;
    E = 73;
else
    sigma = 235;
    ro = 7850;
    E = 200;
end

t = input('input t: ');
a = input('input a: ');
b = input('input b: ');

while true
    c = input('input c: ');
    if c >= 30
        break;
    end
    disp('try again');
end

while true
    R1 = input('input R1: ');
    if R1 >= 3*t
        break;
    end
    disp('try again');
end

while true
    R2 = input('input R2: ');
    if R2 >= 3*t
        break;
    end
    disp('try again');
end

% محاسبات اولیه
H = 2*t + R1 + R2 + a;
Yc = (b*t*t/2 + c*t*(1.5*t+R1+R2) + 2*a*t*(t+R2+a/2) + ...
      pi*R1*t/2*(t+R2+a+2*R1/pi) + pi*R2*t/2*(t+2*R2/pi)) / ...
     (b*t + c*t + 2*a*t + pi*t*0.5*(R1+R2));

I = 1/12*(b*t^3 + c*t^3 + 2*t*a^3) + t*(2/pi)*(R1^3 + R2^3) + ...
    b*t*(t/2-Yc)^2 + c*t*(3*t/2+R1+R2+a-Yc)^2 + 2*a*t*(t+R2+a/2-Yc)^2 + ...
    pi*t/2*(R1*(t+a+R2+2*R1/pi-Yc)^2 + R2*(t+2*R2/pi-Yc)^2);

if H-Yc < Yc
    C = Yc;
else
    C = H-Yc;
end

A = t*(2*b + 2*a - c + pi*(R1+R2+t));
F = C/I;

fprintf('A = %.2f\n', A);
fprintf('F = %.6f\n', F);

Q = 4*F*1e3*(500 + 1.8e-6*ro*A);
Y = (2000/48 + 100*ro*A*(0.0002 + 5/384))/(E*I);

fprintf('Q = %.2f MPa\n', Q);
fprintf('Y = %.6e m\n', Y);

% تابع هدف: مینیمم کردن مساحت
objective = @(vars) evaluate_objective(vars, ro);

% قیود نامساوی
nonlcon = @(vars) nonlinear_constraints(vars, sigma, ro, E);

% مقدار اولیه
x0 = [a, b, c, R1, R2, t];

% محدودیت‌های متغیرها
lb = [5, 5, 30, 3, 3, 1];
ub = [1000, 1000, 1000, 500, 500, 10];

% تنظیمات بهینه‌سازی
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 10000, ...
    'Algorithm', 'sqp', ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

% حل مسئله بهینه‌سازی
[x_opt, fval, exitflag, output] = fmincon(objective, x0, [], [], [], [], ...
    lb, ub, nonlcon, options);

% نمایش نتایج
if exitflag > 0
    a_opt = x_opt(1);
    b_opt = x_opt(2);
    c_opt = x_opt(3);
    R1_opt = x_opt(4);
    R2_opt = x_opt(5);
    t_opt = x_opt(6);
    
    [A_opt, Q_opt, Y_opt] = evaluate(x_opt, sigma, ro, E);
    
    fprintf('\n--- طرح بهینه (حداقل مساحت) ---\n');
    fprintf('a = %.2f, b = %.2f, c = %.2f\n', a_opt, b_opt, c_opt);
    fprintf('R1 = %.2f, R2 = %.2f, t = %.2f\n', R1_opt, R2_opt, t_opt);
    fprintf('مساحت = %.2f mm²\n', A_opt);
    fprintf('Q = %.2f MPa  (باید ≤ %.2f باشد)\n', Q_opt, sigma);
    fprintf('Y = %.6e m  (باید ≤ 5e-6 باشد)\n', Y_opt);
else
    fprintf('بهینه‌سازی موفق نبود. exitflag = %d\n', exitflag);
    disp(output.message);
end

%% توابع کمکی

function A = evaluate_objective(vars, ro)
    t = vars(6);
    a = vars(1);
    b = vars(2);
    c = vars(3);
    R1 = vars(4);
    R2 = vars(5);
    A = t*(2*b + 2*a - c + pi*(R1+R2+t));
end

function [A, Q, Y] = evaluate(vars, sigma, ro, E)
    a = vars(1);
    b = vars(2);
    c = vars(3);
    R1 = vars(4);
    R2 = vars(5);
    t = vars(6);
    
    H = 2*t + R1 + R2 + a;
    
    Yc = (b*t*(t/2) + (pi/2)*R2*t*(t+2*R2/pi) + 2*a*t*(t+R2+a/2) + ...
          (pi/2)*R1*t*(t+R2+a+2*R1/pi) + c*t*(t+R2+a+R1+t/2)) / ...
         (b*t + (pi/2)*R2*t + 2*a*t + (pi/2)*R1*t + c*t);
    
    I = b*t*(t/2-Yc)^2 + ...
        ((2/pi)*R2^3*t + (pi/2)*R2*t*(t+2*R2/pi-Yc)^2) + ...
        2*a*t*(t+R2+a/2-Yc)^2 + ...
        ((2/pi)*R1^3*t + (pi/2)*R1*t*(t+R2+a+2*R1/pi-Yc)^2) + ...
        c*t*(t+R2+a+R1+t/2-Yc)^2;
    
    if H-Yc < Yc
        C = Yc;
    else
        C = H-Yc;
    end
    
    A = t*(2*b + 2*a - c + pi*(R1+R2+t));
    F = C/I;
    Q = 4*F*1e3*(600 + 1.8e-6*ro*A);
    Y = (2000/48 + 100*ro*A*(0.0002 + 5/384))/(E*I);
end

function [c_ineq, c_eq] = nonlinear_constraints(vars, sigma, ro, E)
    a = vars(1);
    b = vars(2);
    c = vars(3);
    R1 = vars(4);
    R2 = vars(5);
    t = vars(6);
    
    [A, Q, Y] = evaluate(vars, sigma, ro, E);
    
    % قیود نامساوی (باید ≤ 0 باشند)
    c_ineq = [
        Q - sigma;              % تنش باید کمتر از حد مجاز باشد
        Y - 5e-6;               % خیز باید کمتر از 5 میلی‌متر باشد
        30 - c;                 % c >= 30
        3*t - R1;               % R1 >= 3*t
        3*t - R2;               % R2 >= 3*t
        c - (b + 2*R1);         % c <= b + 2*R1
    ];
    
    % قیود مساوی (ندارید)
    c_eq = [];
end