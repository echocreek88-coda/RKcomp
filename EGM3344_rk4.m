function [t, y] = rk4(tspan, y0, dy_dt,n)
% RK4    Fourth-order Runge-Kutta for an IVP.
% Input:
%   du_dt   defines f in y'(t) = f(t, u) 
%   tspan   endpoints of time interval (2-vector)
%   y0      initial value (m-vector)
%   n       number of time steps (integer)
% Output:
%   t       selected nodes (vector, length n+1)
%   y       solution values (array, n+1 by m)
    
    if nargin < 3 %default ODE builder
        dy_dt = 0;
        n = tspan(2)*100;
    elseif nargin < 4
        n = tspan(2)*100; % default steps
    end
    
    % Define the time discretization.
    a = tspan(1);  b = tspan(2);
    h = (b - a) / n;
    t = a + (0:n)' * h;

    % Initialize solution array. 
    y = zeros(length(y0), n+1);
    y(:, 1) = y0(:);

    % Time stepping.
    for i = 1:n
      k1 = h * dy_dt( t(i),       y(:, i)       );
      k2 = h * dy_dt( t(i) + h/2, y(:, i) + k1/2);
      k3 = h * dy_dt( t(i) + h/2, y(:, i) + k2/2);
      k4 = h * dy_dt( t(i) + h,   y(:, i) + k3  );
      y(:, i+1) = y(:, i) + (k1 + 2*(k2 + k3) + k4) / 6;
    end

    y = y.';    % conform to MATLAB output convention
end

f = @(t,y) -y.^2;

[t,y]=rk4([0 30], 10, f);