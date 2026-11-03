% انتخاب جنس
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

% ورودی‌های اولیه
t = input('input t: ');
a = input('input a: ');
b = input('input b: ');

% ورودی teta با شرط
while true
    teta = input('input teta: ');
    radian = deg2rad(teta);
    s = sin(radian);
    c = cos(radian);
    if b + 2*a*s > 30 && teta > -90 && teta <= 90
        break;
    end
    disp('try again');
end

c1 = input('input c1: ');

% محاسبات اولیه
radian = deg2rad(teta);
s = sin(radian);
c = cos(radian);

H = t + a*c;
Yc = (c1*t + c*a^2 + (a*c + 0.5*t)*b) / (2*c1 + 2*a + b);
I = (0.5*t - Yc)^2 * (2*c1*t) + (0.5*a*c - Yc)^2 * (2*a*t) + (a*c - Yc)^2 * (b*t) + ...
    (2*c1*t^3 + b*t^3 + 2*(a*s*(t*s)^3) + 2*t*c*(a*c)^3 + b*t^3) / 12;

if H - Yc < Yc
    C = Yc;
else
    C = H - Yc;
end

A = t * (2*c1 + 2*a + b);
F = C / I;

fprintf('A = %.2f\n', A);
fprintf('I = %.2f\n', I);
fprintf('F = %.2f\n', F);

Q = 4 * F * 10^3 * (500 + 1.8 * 10^(-6) * ro * A);
Y = 10^6 * (2000/48 + 100*ro*A*(0.0002 + 5/384)) / (E*I);

fprintf('Q = %.2f\n', Q);
fprintf('Y = %.2f\n', Y);

% تابع هدف
objective = @(vars) objectiveFunc(vars, sigma, ro, E);

% قیود
nonlcon = @(vars) constraintsFunc(vars, sigma, ro, E);

% نقطه شروع
x0 = [a, b, c1, teta, t];

% محدودیت‌های متغیرها
lb = [5, 5, 5, -90, 1];
ub = [10000, 10000, 10000, 90, 10];

% تنظیمات بهینه‌سازی
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 500, ...
    'Algorithm', 'sqp');

% بهینه‌سازی
[x_opt, A_opt, exitflag, output] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);

% نمایش نتایج
if exitflag > 0
    a_opt = x_opt(1);
    b_opt = x_opt(2);
    c1_opt = x_opt(3);
    teta_opt = x_opt(4);
    t_opt = x_opt(5);
    
    [~, Q_opt, Y_opt] = evaluate(x_opt, sigma, ro, E);
    
    fprintf('\n--- طرح بهینه (حداقل مساحت) ---\n');
    fprintf('a = %.2f, b = %.2f, c1 = %.2f\n', a_opt, b_opt, c1_opt);
    fprintf('teta = %.2f, t = %.2f\n', teta_opt, t_opt);
    fprintf('مساحت = %.2f mm²\n', A_opt);
    fprintf('Q = %.2f MPa  (باید ≤ %.2f باشد)\n', Q_opt, sigma);
    fprintf('Y = %.2f mm  (باید ≤ 5 باشد)\n', Y_opt);
else
    fprintf('بهینه‌سازی موفق نبود. exitflag = %d\n', exitflag);
    disp(output.message);
end

%% توابع کمکی

function [A, Q, Y] = evaluate(vars, sigma, ro, E)
    a = vars(1);
    b = vars(2);
    c1 = vars(3);
    teta = vars(4);
    t = vars(5);
    
    radian = deg2rad(teta);
    s = sin(radian);
    c = cos(radian);
    
    H = t + a*c;
    Yc = (c1*t + c*a^2 + (a*c + 0.5*t)*b) / (2*c1 + 2*a + b);
    I = (0.5*t - Yc)^2 * (2*c1*t) + (0.5*a*c - Yc)^2 * (2*a*t) + (a*c - Yc)^2 * (b*t) + ...
        (2*c1*t^3 + b*t^3 + 2*(a*s*(t*s)^3) + 2*t*c*(a*c)^3 + b*t^3) / 12;
    
    if H - Yc < Yc
        C = Yc;
    else
        C = H - Yc;
    end
    
    A = t * (2*c1 + 2*a + b);
    F = C / I;
    Q = 4 * F * 10^3 * (500 + 1.8 * 10^(-6) * ro * A);
    Y = 10^6 * (2000/48 + 100*ro*A*(0.0002 + 5/384)) / (E*I);
end

function f = objectiveFunc(vars, sigma, ro, E)
    [A, ~, ~] = evaluate(vars, sigma, ro, E);
    f = A;
end

function [c, ceq] = constraintsFunc(vars, sigma, ro, E)
    a = vars(1);
    b = vars(2);
    teta = vars(4);
    
    radian = deg2rad(teta);
    s = sin(radian);
    
    [~, Q, Y] = evaluate(vars, sigma, ro, E);
    
    % قیود نامساوی (باید ≤ 0)
    c = [Q - sigma;      % Q باید ≤ sigma
         Y - 5;          % Y باید ≤ 5
         30 - (b + 2*a*s)]; % b + 2*a*s باید > 30
    
    % قیود مساوی (ندارد)
    ceq = [];
end