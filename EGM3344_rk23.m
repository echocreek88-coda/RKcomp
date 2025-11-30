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


function [t, r1,r2,r3] = rk23(tspan, y0, dy_dt, tol)
% RK23   Adaptive IVP solver based on embedded RK formulas.
% Input:
%   du_dt   defines f in u'(t) = f(t, u) 
%   tspan   endpoints of time interval (2-vector)
%   u0      initial value (m-vector)
%   tol     global error target (positive scalar)
% Output:
%   t       selected nodes (vector, length n+1)
%   u       solution values (array, n+1 by m)

    if nargin < 3 %default ODE builder
        dy_dt = @threeBodyODE;
        tol = 1e-6;
    elseif nargin < 4
        tol = 1e-6;  % default tolerance
    end


% Initialize for the first time step.
    i = 1;
    t = tspan(1);
    y(:, 1) = y0(:);
    h = 0.5 * tol^(1/3);
    s1 = dy_dt(t(1), y(:, 1));

    % Time stepping.
    while t(i) < tspan(2)
        % Detect underflow of the step size.
        if t(i) + h == t(i)
            warning('Stepsize too small near t=%.6g.',t(i))
            break  % quit time stepping loop
        end
        
        % New RK stages.
        s2 = dy_dt(t(i) + h/2,   y(:, i) + (h/2)   * s1);
        s3 = dy_dt(t(i) + 3*h/4, y(:, i) + (3*h/4) * s2);
        unew2 = y(:, i) + h * (2*s1 + 3*s2 + 4*s3) / 9;    % 2rd order solution
        s4 = dy_dt(t(i) + h, unew2);
        err = h * (-5*s1/72 + s2/12 + s3/9 - s4/8);        % 2nd/3rd order difference
        E = norm(err, Inf);                                % error estimate
        maxerr = tol * (1 + norm(y(:, i), Inf));           % relative/absolute blend
        
        % Accept the proposed step? 
        if E < maxerr     % yes 
            t = [t; t(i) + h];
            y(:, i+1) = unew2;
            i = i+1;
            s1 = s4;      % use FSAL property
        end
        
        % Adjust step size. 
        q = 0.8 * (maxerr / E) ^ (1/3);    % optimal step factor
        q = min(q, 4);                     % limit stepsize growth
        h = min(q*h, tspan(2) - t(i));     % don't step past the end
    end

    r1 = y(1:3, :)';
    r2 = y(7:9,:)';
    r3 = y(13:15,:)';    % conform to MATLAB output convention
end
%{
r1 = [2,0,0];
v1 = [0,1,0];
r2 = [0,1,0];
v2 = [0,0,0];
r3 = [0,0,1];
v3 = [0,0,0];
%}
y0 = [r1, v1, r2, v2, r3, v3];  % pass in vals as row vectors
tspan = [0 1e3]; %initial and final time value

[t,r1,r2,r3]=rk23(tspan,y0 );