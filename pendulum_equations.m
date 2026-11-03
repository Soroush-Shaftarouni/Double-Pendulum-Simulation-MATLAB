function dz = pendulum_equations(t, z, p)
    teta1 = z(1); omega1 = z(2);
    teta2 = z(3); omega2 = z(4);
    
    delta = teta1 - teta2;
    
    M11 = (p.m1 + p.m2) * p.L1;
    M12 = p.m2 * p.L2 * cos(delta);
    M21 = p.m2 * p.L1 * cos(delta);
    M22 = p.m2 * p.L2;
    
    M = [M11, M12; M21, M22];
    
    B1 = -p.m2 * p.L2 * omega2^2 * sin(delta) - (p.m1 + p.m2) * p.g * sin(teta1);
    B2 =  p.m2 * p.L1 * omega1^2 * sin(delta) - p.m2 * p.g * sin(teta2);
    
    B = [B1; B2];
    
    accel = M \ B;
    
    dz = zeros(4,1);
    dz(1) = omega1;      
    dz(2) = accel(1);    
    dz(3) = omega2;      
    dz(4) = accel(2);    
end

