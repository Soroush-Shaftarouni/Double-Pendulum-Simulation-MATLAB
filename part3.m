clc; clear; close all;

m1 = 1; m2 = 1; g = 9.81;

F = @(x) [(m1+m2)*g*sin(x(1)); m2*g*sin(x(2))];

J = @(x) [(m1+m2)*g*cos(x(1)), 0; 
          0, m2*g*cos(x(2))];

x0 = [0.5; 0.5];

x_nr = x0;
tol = 1e-6;
max_iter = 100;
iter_nr = 0;

for i = 1:max_iter
    Fx = F(x_nr);
    if norm(Fx) < tol
        break;
    end
    Jx = J(x_nr);
    dx = -Jx \ Fx;
    x_nr = x_nr + dx;
    iter_nr = iter_nr + 1;
end

error_nr = norm(F(x_nr));
fprintf('--- نیوتن-رافسون ---\n');
fprintf('ریشه ها (teta1, teta2): [%.6f, %.6f]\n', x_nr(1), x_nr(2));
fprintf('تعداد تکرار: %d\n', iter_nr);
fprintf('خطای نهایی: %e\n\n', error_nr);

options = optimoptions('fsolve', 'Display', 'none');
[x_fs, fval_fs, exitflag, output] = fsolve(F, x0, options);

error_fs = norm(fval_fs);
fprintf('--- fsolve متلب ---\n');
fprintf('ریشه ها (teta1, teta2): [%.6f, %.6f]\n', x_fs(1), x_fs(2));
fprintf('تعداد تکرار: %d\n', output.iterations);
fprintf('خطای نهایی: %e\n', error_fs);