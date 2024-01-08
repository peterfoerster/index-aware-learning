plot_solution    = 0; % 0, 1
plot_consistency = 1; % 0, 1

variables = ['i_L_1'; 'i_L_2'; 'v_1'; 'v_2'; 'v_3'; 'v_4'; 'i_V'];
tol       = [5e-3]; % 1e-2, 5e-3, 1e-3

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

L_1  = 27.46e-6;
L_12 = 27.57e-6;
L_2  = 27.75e-6;

v_s = @(t) 12*sin(2*pi*50*t);

n_v = 4;
n_L = 2;
n_V = 1;

A_R = [[0, 0, 0, 0, 0]; [1, 0, -1, 0, 0]; [-1, -1, 0, 0, 1]; [0, 0, 1, 1, -1]];
A_L = [[1 0]; [0 1]; [0, 0]; [0, 0]];
A_V = [[1]; [0]; [0]; [0]];

G = @(x, T) [[g_D(x(2)-x(3), T(1)), 0, 0, 0, 0]; [0, g_D(-x(3), T(2)), 0, 0, 0]; [0, 0, g_D(x(4)-x(2), T(3)), 0, 0]; [0, 0, 0, g_D(x(4), T(4)), 0]; [0, 0, 0, 0, 1/R]];
L = [[L_1 L_12]; [L_12 L_2]];

M = [[zeros(n_v,n_v), zeros(n_v,n_L), zeros(n_v,n_V)]; [zeros(n_L,n_v), L, zeros(n_L,n_V)]; [zeros(n_V,n_v), zeros(n_V,n_L), zeros(n_V,n_V)];];
K = @(x, T) [[A_R*G(x, T)*A_R.', A_L, A_V]; [-A_L.', zeros(n_L,n_L), zeros(n_L,n_V)]; [-A_V.', zeros(n_V,n_L), zeros(n_V,n_V)]];
f = @(t) [zeros(n_v,1); zeros(n_L,1); v_s(t)];

x_0 = zeros(n_v+n_L+n_V,1);

% evaluation point
T = [65; 85; 85; 65] + 273.15;

if (plot_solution)
    for iv = 1:length(variables)
        variable = variables(iv,:);
        filename = ['fwr_Deltat=' num2str(Deltat) '_T_1=' num2str(T(1)) '_T_2=' num2str(T(2)) '_T_3=' num2str(T(3)) '_T_4=' num2str(T(4)) '.dat'];
        if (~exist(filename))
            % [x_d, t_d] = trapezoidal_rule(t_0, t_f, Deltat, M, @(x) K(x, T), @(t) f(t), x_0);
            [x_d, t_d] = implicit_euler(t_0, t_f, Deltat, M, @(x) K(x, T), @(t) -f(t), x_0);
            write_dat (filename, t_d, x_d);
        else
            data = dlmread(filename);
            data = data(2:end,:);
            t_d  = data(:,1)';
            x_d  = data(:,2:end)';
        end

        if (strcmp(variable, 'i_L_1'))
            fprintf('\n%s\n', variable);
            x = interp1(t_d, x_d(5,:), t);
        elseif (strcmp(variable, 'i_L_2'))
            fprintf('\n%s\n', variable);
            x = interp1(t_d, x_d(6,:), t);
        else
            variable = variable(1:3);
            if (strcmp(variable, 'v_1'))
                fprintf('\n%s\n', variable);
                x = interp1(t_d, x_d(1,:), t);
            elseif (strcmp(variable, 'v_2'))
                fprintf('\n%s\n', variable);
                x = interp1(t_d, x_d(2,:), t);
            elseif (strcmp(variable, 'v_3'))
                fprintf('\n%s\n', variable);
                x = interp1(t_d, x_d(3,:), t);
            elseif (strcmp(variable, 'v_4'))
                fprintf('\n%s\n', variable);
                x = interp1(t_d, x_d(4,:), t);
            elseif (strcmp(variable, 'i_V'))
                fprintf('\n%s\n', variable);
                x = interp1(t_d, x_d(7,:), t);
            end
        end
        x_p = x';

        % all sampling points
        % T_1 = T_2 = T_3 = T_4 = [20:10:100] + 273.15;
        T_1 = T_2 = [20:10:100] + 273.15;

        N_T_1 = length(T_1);
        N_T_2 = length(T_2);
        % N_T_3 = length(T_3);
        % N_T_4 = length(T_4);
        X     = NaN(N_t*N_T_1*N_T_2, 3);
        y     = NaN(N_t*N_T_1*N_T_2, 1);
        for i_2 = 1:N_T_2
            for i_1 = 1:N_T_1
                filename = ['fwr_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1(i_1)) '_T_2=' num2str(T_2(i_2)) '_T_3=' num2str(T_2(i_2)) '_T_4=' num2str(T_1(i_1)) '.dat'];
                data     = dlmread(filename);
                data     = data(2:end,:);
                t_d      = data(:,1)';
                x_d      = data(:,2:end)';

                i_t                        = (i_1 - 1) + (i_2 - 1)*N_T_1;
                X(i_t*N_t+1:(i_t+1)*N_t,:) = [t', T_1(i_1)*ones(N_t,1), T_2(i_2)*ones(N_t,1)];


                if (strcmp(variable, 'i_L_1'))
                    x = interp1(t_d, x_d(5,:), t);
                elseif (strcmp(variable, 'i_L_2'))
                    x = interp1(t_d, x_d(6,:), t);
                else
                    variable = variable(1:3);
                    if (strcmp(variable, 'v_1'))
                        x = interp1(t_d, x_d(1,:), t);
                    elseif (strcmp(variable, 'v_2'))
                        x = interp1(t_d, x_d(2,:), t);
                    elseif (strcmp(variable, 'v_3'))
                        x = interp1(t_d, x_d(3,:), t);
                    elseif (strcmp(variable, 'v_4'))
                        x = interp1(t_d, x_d(4,:), t);
                    elseif (strcmp(variable, 'i_V'))
                        x = interp1(t_d, x_d(7,:), t);
                    end
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

            filename = ['fwr_I_' variable '_tol=' num2str(tol(i_tol)) '.dat'];
            I        = dlmread(filename);

            filename = ['fwr_param_' variable '_tol=' num2str(tol(i_tol)) '.dat'];
            GP.param = dlmread(filename);

            filename            = ['fwr_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '.dat'];
            GP.lognoisevariance = dlmread(filename);

            tic;
            y_p = stk_predict(GP, X(I,:), y(I), X_p);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p(i_tol) = norm(y_p.mean - x_p) / norm(x_p)

            filename = ['fwr_solution_' variable '_tol=' num2str(tol(i_tol)) '.dat'];
            write_dat(filename, t, [x_p'; y_p.mean']);
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
