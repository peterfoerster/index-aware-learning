plot_convergence = 1; % 0, 1
plot_solution    = 0; % 0, 1
plot_consistency = 0; % 0, 1

variables = ['i_L_2'; 'phi_1'; 'phi_2'; 'phi_3'; 'phi_4'; 'i_L_1'];
tol       = [1e-3]; % 1e-2, 5e-3, 1e-3

t_0       = 0;
t_f       = 5e-2;
Deltat    = 1e-5;
N_t       = 1000;
t         = linspace(t_0, t_f, N_t);
epsilon_T = 1e-6;

% T in K
q   = 1.602e-19;
k   = 1.38e-23;
qk  = q/k;
I_s = @(T) 1e-6 * (T./300).^3 .* exp( (T./300 - 1) .* 1.12*qk./T );
g_D = @(v_D, T) I_s(T) .* qk./T .* exp(qk./T .* v_D);
R   = 50;

% lullo2017: (4) and (7)
L_m  = @(i_L) 3.5e-3 - 1.2e-4*i_L - 1.3e-4*i_L.^2 + 3e-5*i_L.^3 - 1.9e-6*i_L.^4;
L_1  = @(i_L_1) L_m(i_L_1);
L_2  = @(i_L_2) 0.1*L_m(i_L_2);
% coupling coefficient k=0.9
L_12 = @(i_L_1, i_L_2) sqrt(0.9*L_m(i_L_1)*0.1*L_m(i_L_2));

i_s  = @(t) 0.1*sin(2*pi*50*t);
i_sp = @(t) 0.1*2*pi*50*cos(2*pi*50*t);

n_phi = 4;
n_L   = 2;
n_V   = 0;

A_R = [[0, 0, 0, 0, 0]; [1, 0, -1, 0, 0]; [-1, -1, 0, 0, 1]; [0, 0, 1, 1, -1]];
A_L = [[1 0]; [0 -1]; [0, 0]; [0, 0]];
A_V = [];
A_I = [-1; 0; 0; 0];

G = @(x, T) [[g_D(x(2)-x(3), T(1)), 0, 0, 0, 0]; [0, g_D(-x(3), T(2)), 0, 0, 0]; [0, 0, g_D(x(4)-x(2), T(3)), 0, 0]; [0, 0, 0, g_D(x(4), T(4)), 0]; [0, 0, 0, 0, 1/R]];
L = @(x) [[L_1(x(5)) L_12(x(5), x(6))]; [L_12(x(5), x(6)) L_2(x(6))]];

M = @(x) [[zeros(n_phi, n_phi), zeros(n_phi, n_L), zeros(n_phi, n_V)]; [zeros(n_L, n_phi), L(x), zeros(n_L, n_V)]; [zeros(n_V, n_phi), zeros(n_V, n_L), zeros(n_V, n_V)]];
K = @(x, T) [[A_R*G(x, T)*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]];
f = @(t) [A_I*i_s(t); zeros(n_L, 1); zeros(n_V, 1)];

x_0 = zeros(n_phi+n_L+n_V,1);

% evaluation point
T = [65; 85; 85; 65] + 273.15;

if (plot_convergence)
    tol = [1e-2, 5e-3, 1e-3];

    for iv = 1:length(variables)
        variable = variables(iv,:);
        N_s      = NaN(length(tol), 1);
        for it = 1:length(tol)
            filename = ['fwr2_I_' variable '_tol=' num2str(tol(it)) '.csv'];
            I        = csvread(filename);
            N_s(it)  = length(I);
        end

        hold on
        semilogy(N_s, tol)
        hold off
        pause

        filename = ['fwr2_convergence_' variable '.csv'];
        csvwrite(filename, [N_s, tol']);
    end
end

if (plot_solution)
    % reference solution
    for iv = 1:length(variables)
        variable = variables(iv,:);
        filename = ['fwr2_Deltat=' num2str(Deltat) '_T_1=' num2str(T(1)) '_T_2=' num2str(T(2)) '_T_3=' num2str(T(3)) '_T_4=' num2str(T(4)) '.csv'];
        if (~exist(filename))
            [x_d, t_d] = implicit_euler(t_0-2*Deltat, t_f, Deltat, @(x) M(x), @(x) K(x, T), @(t) -f(t), x_0);
            csvwrite(filename, [t_d; x_d]');
        end
        data = csvread(filename);
        % exclude first two timesteps
        data = data(3:end,:);
        t_d  = data(:,1)';
        x_d  = data(:,2:end)';

        if (strcmp(variable, 'i_L_2'))
            x = interp1(t_d, x_d(6,:), t);
        elseif (strcmp(variable, 'phi_1'))
            x = interp1(t_d, x_d(1,:), t);
        elseif (strcmp(variable, 'phi_2'))
            x = interp1(t_d, x_d(2,:), t);
        elseif (strcmp(variable, 'phi_3'))
            x = interp1(t_d, x_d(3,:), t);
        elseif (strcmp(variable, 'phi_4'))
            x = interp1(t_d, x_d(4,:), t);
        elseif (strcmp(variable, 'i_L_1'))
            x = interp1(t_d, x_d(5,:), t);
        end
        x_p = x';

        hold on
        plot(t, x)
        hold off
        pause

        % all sampling points
        T_1 = T_2 = [20:10:100] + 273.15;

        N_T_1 = length(T_1);
        N_T_2 = length(T_2);
        X     = NaN(N_t*N_T_1*N_T_2, 3);
        y     = NaN(N_t*N_T_1*N_T_2, 1);
        for i_2 = 1:N_T_2
            for i_1 = 1:N_T_1
                filename = ['fwr2_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1(i_1)) '_T_2=' num2str(T_2(i_2)) '_T_3=' num2str(T_2(i_2)) '_T_4=' num2str(T_1(i_1)) '.csv'];
                data     = csvread(filename);
                % exclude first two timesteps
                data     = data(3:end,:);
                t_d      = data(:,1)';
                x_d      = data(:,2:end)';

                i_t                        = (i_1 - 1) + (i_2 - 1)*N_T_1;
                X(i_t*N_t+1:(i_t+1)*N_t,:) = [t', T_1(i_1)*ones(N_t,1), T_2(i_2)*ones(N_t,1)];

                if (strcmp(variable, 'i_L_2'))
                    x = interp1(t_d, x_d(6,:), t);
                elseif (strcmp(variable, 'phi_1'))
                    x = interp1(t_d, x_d(1,:), t);
                elseif (strcmp(variable, 'phi_2'))
                    x = interp1(t_d, x_d(2,:), t);
                elseif (strcmp(variable, 'phi_3'))
                    x = interp1(t_d, x_d(3,:), t);
                elseif (strcmp(variable, 'phi_4'))
                    x = interp1(t_d, x_d(4,:), t);
                elseif (strcmp(variable, 'i_L_1'))
                    x = interp1(t_d, x_d(5,:), t);
                end
                y(i_t*N_t+1:(i_t+1)*N_t) = x';

                % check if solution is parameter dependent
                if (i_1 == i_2 == 1)
                    xh = x;
                else
                    delta_T = norm(x - xh) / norm(xh);
                    xh = x;
                end
            end
        end

        % only learn time dependence if solution is not parameter dependent
        if (delta_T < epsilon_T)
            fprintf('\n%s does not depend on parameters\n', variable);

            % only consider one set of temperatures
            X = X(1:N_t,1);
            y = y(1:N_t);

            X_p = t';
            GP  = stk_model('stk_gausscov_aniso', 1);
        else
            X_p = [t', T(1)*ones(length(t), 1), T(2)*ones(length(t), 1)];
            GP  = stk_model('stk_gausscov_aniso', 3);
        end

        e_p = NaN(length(tol), 1);
        for i_tol = 1:length(tol)
            fprintf('\ntol: %g\n', tol(i_tol));

            filename = ['fwr2_I_' variable '_tol=' num2str(tol(i_tol)) '.csv'];
            I        = csvread(filename);

            filename = ['fwr2_param_' variable '_tol=' num2str(tol(i_tol)) '.csv'];
            GP.param = csvread(filename);

            filename            = ['fwr2_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '.csv'];
            GP.lognoisevariance = csvread(filename);

            tic;
            y_p = stk_predict(GP, X(I,:), y(I), X_p);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p(i_tol) = norm(y_p.mean - x_p) / norm(x_p)

            filename = ['fwr2_solution_' variable '_tol=' num2str(tol(i_tol)) '.csv'];
            csvwrite(filename, [t; x_p'; y_p.mean']');
        end
    end
end

if (plot_consistency)
    % first basis functions
    Q      = W      = zeros(7,5);
    Q(1,1) = W(1,1) = 1;
    Q(2,2) = W(2,2) = 1;
    Q(3,3) = W(3,3) = 1;
    Q(4,4) = W(4,4) = 1;
    Q(7,5) = W(7,5) = 1;
    P      = V      = zeros(7,2);
    P(5,1) = V(5,1) = 1;
    P(6,2) = V(6,2) = 1;

    % first stage
    Mt   = V.' * M * P;
    Kt_P = @(x,T) V.' * K(x, T) * P;
    Kt_Q = @(x,T) V.' * K(x, T) * Q;
    ft   = @(t) V.' * f(t);
    Kb_P = @(x,T) W.' * K(x, T) * P;
    Kb_Q = @(x,T) W.' * K(x, T) * Q;
    fb   = @(t) W.' * f(t);

    % algebraic function
    g = @(xt, xb, t) Kb_P(P*xt + Q*xb, T) * xt + Kb_Q(P*xt + Q*xb, T) * xb + fb(t);
    % options.TolX   = 1e-16;
    % options.TolFun = 1e-16;
    % xb = fsolve(@(xb) g(xt, xb, t), xb, options);

    % [t, i_L_1, i_L_1_p]
    datname = ['fwr_solution_' 'i_L_1' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    t       = data(2:end,1);
    i_L_1   = data(2:end,2);
    i_L_1_p = data(2:end,3);

    datname = ['fwr_solution_' 'i_L_2' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    i_L_2   = data(2:end,2);
    i_L_2_p = data(2:end,3);

    datname = ['fwr_solution_' 'v_1' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    v_1     = data(2:end,2);
    v_1_p   = data(2:end,3);

    datname = ['fwr_solution_' 'v_2' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    v_2     = data(2:end,2);
    v_2_p   = data(2:end,3);

    datname = ['fwr_solution_' 'v_3' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    v_3     = data(2:end,2);
    v_3_p   = data(2:end,3);

    datname = ['fwr_solution_' 'v_4' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    v_4     = data(2:end,2);
    v_4_p   = data(2:end,3);

    datname = ['fwr_solution_' 'i_V' '_tol=' num2str(tol) '.dat'];
    data    = dlmread(datname);
    i_V     = data(2:end,2);
    i_V_p   = data(2:end,3);

    xt   = [i_L_1'; i_L_2'];
    xt_p = [i_L_1_p'; i_L_2_p'];
    xb   = [v_1'; v_2'; v_3'; v_4'; i_V'];
    xb_p = [v_1_p'; v_2_p'; v_3_p'; v_4_p'; i_V_p'];
    xb_r = NaN(size(xb_p));
    e_s  = NaN(1,length(t));
    e_c  = NaN(1,length(t));
    e_r  = NaN(1,length(t));
    options.TolX   = 1e-16;
    options.TolFun = 1e-16;
    for it = 1:length(t)
        if (mod(it, 1000) == 0)
            fprintf('\ntimestep no.: %i/%i\n', it, length(t));
        end

        % start with initial condition or previous timestep
        if it == 1
            xb_r(:,it) = fsolve(@(xb) g(xt_p(:,it), xb, t(it)), xb(:,it), options);
        else
            xb_r(:,it) = fsolve(@(xb) g(xt_p(:,it), xb, t(it)), xb_r(:,it-1), options);
        end

        e_s(it) = norm(g(xt(:,it), xb(:,it), t(it)));
        e_c(it) = norm(g(xt_p(:,it), xb_p(:,it), t(it)));
        e_r(it) = norm(g(xt_p(:,it), xb_r(:,it), t(it)));
    end

    max(e_s)
    max(e_c)
    max(e_r)

    semilogy(t, e_s, t, e_c, t, e_r);

    filename = 'fwr_solution.dat';
    write_dat (filename, t', [xt; xt_p; xb; xb_p; xb_r]);

    filename = 'fwr_consistency.dat';
    write_dat (filename, t', [e_s; e_c; e_r]);
end
