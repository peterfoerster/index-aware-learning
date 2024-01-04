%   [x, t] = implicit_euler (T_0, T, Deltat, M, K, f, x_0)
%
% Computes the numerical solution of an initial value problem of the form
%      M(x) x' + K(x) x = f(t) on [T_0,T]
%               x(0) = x_0
% using the implicit Euler method.
%
% INPUT:
%       - T_0, T, Deltat: time interval [T_0,T] and timestep
%       - M: @(x) nonlinear matrix
%       - K: @(x) nonlinear matrix
%       - f: @(t) time dependent rhs
%       - x_0: initial condition
%
% OUTPUT:
%       - x: [N_x, N_t] solution vector at the discrete points in time
%       - t: [1, N_t] vector containing the discrete points in time

function [x, t] = implicit_euler (T_0, T, Deltat, M, K, f, x_0)
    N_x = length(x_0);
    N_t = round((T - T_0) / Deltat);
    x   = NaN(N_x, N_t+1);
    x(:,1) = x_0;

    t = T_0:Deltat:(T_0 + N_t*Deltat);
    for n=1:N_t
        if (~mod(n, 1000))
            fprintf('\niteration no.: %i/%i\n', n, N_t);
        end

        % M(x_npo) (x_npo - x_n) + Deltat K(x_npo) x_npo = Deltat f(t_npo)
        g = @(x_npo) M(x_npo)*(x_npo - x(:,n)) + Deltat*K(x_npo)*x_npo - Deltat*f(t(n+1));

        options.TolX   = 1e-16;
        options.TolFun = 1e-16;
        x(:,n+1)       = fsolve(g, x(:,n), options);
    end
end
