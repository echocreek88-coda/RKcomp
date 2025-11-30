function dydt = threeBodyODE(t,y)
% threeBodyODE_simple - explicit 3-body ODE system for use with myODE45
%
% y = [x1 y1 z1 vx1 vy1 vz1 x2 y2 z2 vx2 vy2 vz2 x3 y3 z3 vx3 vy3 vz3]'

    G = 6.674e-11;           % gravitational constant (can be scaled)
    m1=1; m2 = 1;  m3 = 1; %masses

    % Extract positions and velocities
    x1 = y(1);  y1 = y(2);  z1 = y(3);
    vx1 = y(4); vy1 = y(5); vz1 = y(6);

    x2 = y(7);  y2 = y(8);  z2 = y(9);
    vx2 = y(10); vy2 = y(11); vz2 = y(12);

    x3 = y(13); y3 = y(14); z3 = y(15);
    vx3 = y(16); vy3 = y(17); vz3 = y(18);

    % Compute distance vectors and magnitudes
    r12 = [x2 - x1; y2 - y1; z2 - z1];
    r13 = [x3 - x1; y3 - y1; z3 - z1];
    r23 = [x3 - x2; y3 - y2; z3 - z2];

    d12 = norm(r12+eps^2*zeros(3,1));
    d13 = norm(r13+eps^2*zeros(3,1));
    d23 = norm(r23+eps^2*zeros(3,1));

    % Accelerations (Newton’s law of gravitation)
    a1 = G*(m2*r12/d12^3 + m3*r13/d13^3);
    a2 = G*(m1*(-r12)/d12^3 + m3*r23/d23^3);
    a3 = G*(m1*(-r13)/d13^3 + m2*(-r23)/d23^3);

    % Collect derivatives
    dydt = [...
        vx1; vy1; vz1; a1(1); a1(2); a1(3); ...
        vx2; vy2; vz2; a2(1); a2(2); a2(3); ...
        vx3; vy3; vz3; a3(1); a3(2); a3(3) ];
end

function [t, r1,r2,r3] = myODE45(tspan, y0, f, tol)
    % myODE45: simple implementation of Runge–Kutta–Fehlberg 4(5)
    % f: function handle @(t,y)
    % tspan: [t0 tf]
    % y0: initial value (scalar or vector)
    % tol: error tolerance

    
    if nargin < 3 %default ODE builder
        f = @threeBodyODE;
        tol = 1e-6;
    elseif nargin < 4
        tol = 1e-6;  % default tolerance
    end

    t0 = tspan(1);
    tf = tspan(2);
    h = (tf - t0)/100;  % initial step guess

    % Coefficients for Dormand–Prince method
    c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    a = [...
        0               0           0           0           0           0;
        1/5             0           0           0           0           0;
        3/40            9/40        0           0           0           0;
        44/45          -56/15       32/9        0           0           0;
        19372/6561    -25360/2187   64448/6561 -212/729     0           0;
        9017/3168     -355/33       46732/5247  49/176    -5103/18656  0;
        35/384          0           500/1113    125/192   -2187/6784    11/84];
    b4 = a(7,:);  % 4th order
    b5 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]; % 5th order
    t = t0;
    y = y0(:);  % store row-wise
    ti = t0;
    yi = y0(:);

    while t(end) < tf
        if t(end) + h > tf
            h = tf - t(end);
        end

        % Compute RK stages
        k1 = f(ti, yi);
        k2 = f(ti + c(2)*h, yi + h*(a(2,1)*k1));
        k3 = f(ti + c(3)*h, yi + h*(a(3,1)*k1 + a(3,2)*k2));
        k4 = f(ti + c(4)*h, yi + h*(a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3));
        k5 = f(ti + c(5)*h, yi + h*(a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4));
        k6 = f(ti + c(6)*h, yi + h*(a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5));

        y4 = yi + h*(b4(1)*k1 + b4(3)*k3 + b4(4)*k4 + b4(5)*k5 + b4(6)*k6);
        y5 = yi + h*(b5(1)*k1 + b5(3)*k3 + b5(4)*k4 + b5(5)*k5 + b5(6)*k6 + b5(7)*f(ti + h, y4));

        % Estimate error
        err = norm(y5 - y4, inf);
        if err < tol
            % Accept step
            ti = ti + h;
            yi = y5;
            t(end+1,1) = ti;
            y = [y yi];
        end

        % Update step size (safety factor 0.9)
        if err == 0
            s = 2;
        else
            s = 0.9 * (tol/err)^(1/5);
        end
        h = h * min(max(s, 0.2), 5);
    end
    r1 = y(1:3,:)';
    r2 = y(7:9,:)';
    r3 = y(13:15,:)';
end



% Test case
%{
r1 = [1,0,0];
v1 = [0,0,0];
r2 = [0,1,0];
v2 = [0,0,0];
r3 = [0,0,1];
v3 = [0,0,0];

y0 = [r1, v1, r2, v2, r3, v3];  % pass in vals as row vectors
tspan = [0 259200/2]; %initial and final time value
%}

% RK45 call
[t, r1, r2, r3] = myODE45(tspan, y0);
%^^^^^^
% Pass in t as a two element row vector [t initial, t final
% Pass y0 as an 6 element row vector [position body 1, velocity body 1, 
% position body 2, velocity body 2, position body 3, velocity body 3] with
% each position/velocity inputted as row vector [x0,y0,z0]








