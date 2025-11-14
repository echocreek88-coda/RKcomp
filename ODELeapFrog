function [t, r1, r2, r3] = myLeapfrog3Body(tspan, y0, f, tol)
% myLeapfrog3Body: 2nd-order symplectic (leapfrog / velocity Verlet)
% Same input/output structure as myODE45, but fixed step, energy-conserving.
%
% tspan: [t0 tf]
% y0   : [x1 y1 z1 vx1 vy1 vz1 x2 y2 z2 vx2 vy2 vz2 x3 y3 z3 vx3 vy3 vz3]
% f    : unused here (kept for interface compatibility)
% tol  : unused here (kept for interface compatibility)
%
% Output:
%   t  : (N+1)x1 time vector
%   r1 : (N+1)x3 positions of body 1
%   r2 : (N+1)x3 positions of body 2
%   r3 : (N+1)x3 positions of body 3

    if nargin < 3
        f   = @threeBodyODE; %#ok<NASGU>  % not used, but kept for compatibility
        tol = 1e-6;                      %#ok<NASGU>
    elseif nargin < 4
        tol = 1e-6;                      %#ok<NASGU>
    end

    t0 = tspan(1);
    tf = tspan(2);

    % Choose fixed number of steps for leapfrog (tweak if needed)
    N  = 2000;
    h  = (tf - t0)/N;

    % Ensure column vector state
    y0 = y0(:);

    % Unpack positions and velocities into 9x1 blocks
    r = [y0(1:3);  y0(7:9);  y0(13:15)];
    v = [y0(4:6);  y0(10:12); y0(16:18)];

    % Initial acceleration
    a = accelFromR(r);

    % Allocate outputs
    t  = zeros(N+1,1);
    r1 = zeros(N+1,3);
    r2 = zeros(N+1,3);
    r3 = zeros(N+1,3);

    % Store initial state
    t(1)    = t0;
    r1(1,:) = r(1:3).';
    r2(1,:) = r(4:6).';
    r3(1,:) = r(7:9).';

    ti = t0;

    % Leapfrog / velocity Verlet loop
    for n = 1:N
        % Half-kick
        v_half = v + 0.5*h*a;

        % Drift
        r_new = r + h*v_half;

        % New acceleration
        a_new = accelFromR(r_new);

        % Half-kick
        v_new = v_half + 0.5*h*a_new;

        % Advance time
        ti = ti + h;

        % Store positions
        t(n+1)    = ti;
        r1(n+1,:) = r_new(1:3).';
        r2(n+1,:) = r_new(4:6).';
        r3(n+1,:) = r_new(7:9).';

        % Update for next step
        r = r_new;
        v = v_new;
        a = a_new;
    end
end

% -------------------------------------------------------------------------
function a = accelFromR(r)
% accelFromR
% Input r = [x1 y1 z1 x2 y2 z2 x3 y3 z3]' (9x1)
% Output a = [ax1 ay1 az1 ax2 ay2 az2 ax3 ay3 az3]' (9x1)
    G  = 1;
    m1 = 1; m2 = 1; m3 = 1;

    x1 = r(1);  y1 = r(2);  z1 = r(3);
    x2 = r(4);  y2 = r(5);  z2 = r(6);
    x3 = r(7);  y3 = r(8);  z3 = r(9);

    r12 = [x2 - x1; y2 - y1; z2 - z1];
    r13 = [x3 - x1; y3 - y1; z3 - z1];
    r23 = [x3 - x2; y3 - y2; z3 - z2];

    d12 = norm(r12 + eps^2*zeros(3,1));
    d13 = norm(r13 + eps^2*zeros(3,1));
    d23 = norm(r23 + eps^2*zeros(3,1));

    a1 = G*( m2*r12/d12^3 + m3*r13/d13^3 );
    a2 = G*( m1*(-r12)/d12^3 + m3*r23/d23^3 );
    a3 = G*( m1*(-r13)/d13^3 + m2*(-r23)/d23^3 );

    a = [a1; a2; a3];
end
