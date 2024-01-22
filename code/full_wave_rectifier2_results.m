plot_convergence = 0; % 0, 1
plot_solution    = 0; % 0, 1
plot_consistency = 1; % 0, 1

variables = ['i_L_2'; 'v_34'; 'phi_1'; 'phi_2'; 'phi_3'; 'phi_4'; 'i_L_1'];
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
C   = 1e-3;

% lullo2017: (4) and (7)
L_m  = @(i_L) 3.5e-3 - 1.2e-4*i_L - 1.3e-4*i_L.^2 + 3e-5*i_L.^3 - 1.9e-6*i_L.^4;
L_1  = @(i_L_1) L_m(i_L_1);
L_2  = @(i_L_2) 0.1*L_m(i_L_2);
% coupling coefficient k=0.9
L_12 = @(i_L_1, i_L_2) sqrt(0.9*L_m(i_L_1)*0.1*L_m(i_L_2));

i_s  = @(t) 0.1*cos(2*pi*50*t);
i_sp = @(t) -0.1*2*pi*50*sin(2*pi*50*t);

n_phi = 4;
n_L   = 2;
n_V   = 0;

A_R = [[0, 0, 0, 0, 0]; [1, 0, -1, 0, 0]; [-1, -1, 0, 0, 1]; [0, 0, 1, 1, -1]];
A_L = [[1 0]; [0 -1]; [0, 0]; [0, 0]];
A_C = [0; 0; 1; -1];
A_V = [];
A_I = [-1; 0; 0; 0];

G = @(x, T) [[g_D(x(2)-x(3), T(1)), 0, 0, 0, 0]; [0, g_D(-x(3), T(2)), 0, 0, 0]; [0, 0, g_D(x(4)-x(2), T(3)), 0, 0]; [0, 0, 0, g_D(x(4), T(4)), 0]; [0, 0, 0, 0, 1/R]];
L = @(x) [[L_1(x(5)) L_12(x(5), x(6))]; [L_12(x(5), x(6)) L_2(x(6))]];

M = @(x) [[A_C*C*A_C.', zeros(n_phi, n_L), zeros(n_phi, n_V)]; [zeros(n_L, n_phi), L(x), zeros(n_L, n_V)]; [zeros(n_V, n_phi), zeros(n_V, n_L), zeros(n_V, n_V)]];
K = @(x, T) [[A_R*G(x, T)*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]];
f = @(t) [A_I*i_s(t); zeros(n_L, 1); zeros(n_V, 1)];

% first basis functions
Q      = W      = zeros(6,3);
Q(1,1) = W(1,1) = 1;
Q(2,2) = W(2,2) = 1;
Q(3,3) = W(3,3) = 1;
Q(4,3) = W(4,3) = 1;
P      = V      = zeros(6,3);
P(3,1) = V(3,1) = 1;
P(5,2) = V(5,2) = 1;
P(6,3) = V(6,3) = 1;

Mt   = @(x) V.' * M(x)*P;
Kt_P = @(x, T) V.' * K(x, T)*P;
Kt_Q = @(x, T) V.' * K(x, T)*Q;
ft   = @(t) V.' * f(t);
Kb_P = @(x, T) W.' * K(x, T)*P;
Kb_Q = @(x, T) W.' * K(x, T)*Q;
fb   = @(t) W.' * f(t);

% second basis functions
Qb      = Wb      = zeros(3,1);
Qb(1,1) = Wb(1,1) = 1;
Pb      = Vb      = zeros(3,2);
Pb(2,1) = Vb(2,1) = 1;
Pb(3,2) = Vb(3,2) = 1;

Qt      = zeros(3,2);
Qt(1,1) = 1;
Qt(3,2) = 1;
Pt      = zeros(3,1);
Pt(2,1) = 1;

Wt = @(x) [0; -L_2(x(6)); L_12(x(5), x(6))];
Vt = @(x) [[1, 0]; [0, L_2(x(6))]; [0, L_12(x(5), x(6))]];

% evaluation point
T = [45; 75; 75; 45] + 273.15;
% phi_3 - phi_4, i_L_2
xt_Q_0  = [0; 0];
% i_L_1
xt_P_0  = 0;
% phi_1
xb_Q_0  = 0;
% phi_2, phi_4
xb_P_0  = [0; 0];

function [g] = g(xt_Q, xt_P, xb_Q, xb_P, t, T, Deltat, Q, P, Qb, Pb, Wb, Vb, Qt, Pt, Wt, Mt, Kt_P, Kt_Q, ft, Kb_P, Kb_Q, fb)
    % x = P xt + Q xb = P (Pt xt_P + Qt xt_Q) + Q (Pb xb_P + Qb xb_Q)
    x       = P*(Pt*xt_P + Qt*xt_Q) + Q*(Pb*xb_P + Qb*xb_Q);

    % approximate derivative by backward difference (assume matrix constant)
    % xt_P    = -(Wb.' * Kb_P(x, T)*Pt) \ (Wb.' * fb(t));
    xt_P_po = -(Wb.' * Kb_P(x, T)*Pt) \ (Wb.' * fb(t+Deltat));
    xt_Pp   = (xt_P_po - xt_P) / Deltat;

    fb_P   = -(Vb.' * Kb_Q(x, T)*Pb) \ (Vb.' * Kb_P(x, T)*Pt*xt_P + Vb.' * fb(t));
    fh     = Mt(x)*Pt*xt_Pp + Kt_P(x, T)*Pt*xt_P + Kt_Q(x, T)*Pb*fb_P + ft(t);
    Kbxt_Q = -(Vb.' * Kb_Q(x, T)*Pb) \ (Vb.' * Kb_P(x, T)*Qt*xt_Q);

    % solve joint nonlinear problem for all algebraic DOFs
    g_xt_P = Wb.' * Kb_P(x, T)*Pt*xt_P + Wb.' * fb(t);
    g_xb_Q = Wt(x).' * Kt_Q(x, T)*Qb*xb_Q + Wt(x).' * (Kt_P(x, T)*Qt*xt_Q + Kt_Q(x, T)*Pb*Kbxt_Q) + Wt(x).' * fh;
    g_xb_P = Vb.' * Kb_Q(x, T)*Pb*xb_P + Vb.' * Kb_P(x, T)*(Qt*xt_Q + Pt*xt_P) + Vb.' * fb(t);

    g = [g_xt_P; g_xb_Q; g_xb_P];
end

% determine consistent initial conditions
g_i            = @(xt_Q, xt_P, xb_Q, xb_P, t) g(xt_Q, xt_P, xb_Q, xb_P, t, T, Deltat, Q, P, Qb, Pb, Wb, Vb, Qt, Pt, Wt, Mt, Kt_P, Kt_Q, ft, Kb_P, Kb_Q, fb);
options.TolX   = 1e-16;
options.TolFun = 1e-16;
x_0            = fsolve(@(x) g_i(xt_Q_0, x(1), x(2), x(3:4), t_0), [xt_P_0; xb_Q_0; xb_P_0], options);
xt_P_0         = x_0(1);
xb_Q_0         = x_0(2);
xb_P_0         = x_0(3:4);
x_0            = P*(Pt*xt_P_0 + Qt*xt_Q_0) + Q*(Pb*xb_P_0 + Qb*xb_Q_0);

if (plot_convergence)
    tol = [1e-2, 5e-3, 1e-3];

    for iv = 1:length(variables)
        variable = variables(iv,:);
        if (variable(1) == 'v')
            variable = variable(1:4);
        end
        N_s = NaN(length(tol), 1);
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
            [x_d, t_d] = implicit_euler(t_0, t_f, Deltat, @(x) M(x), @(x) K(x, T), @(t) -f(t), x_0);
            csvwrite(filename, [t_d; x_d]');
        end
        data = csvread(filename);
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
        else
            variable = variable(1:4);
            if (strcmp(variable, 'v_34'))
                x = interp1(t_d, x_d(3,:)-x_d(4,:), t);
            end
        end
        x_p = x';

        hold on
        plot(t, x_p)
        hold off
        pause

        % all sampling points
        T_1 = T_2 = [20:10:90] + 273.15;

        N_T_1 = length(T_1);
        N_T_2 = length(T_2);
        X     = NaN(N_t*N_T_1*N_T_2, 3);
        y     = NaN(N_t*N_T_1*N_T_2, 1);
        for i_2 = 1:N_T_2
            for i_1 = 1:N_T_1
                filename = ['fwr2_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1(i_1)) '_T_2=' num2str(T_2(i_2)) '_T_3=' num2str(T_2(i_2)) '_T_4=' num2str(T_1(i_1)) '.csv'];
                data     = csvread(filename);
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
                else
                    variable = variable(1:4);
                    if (strcmp(variable, 'v_34'))
                        x = interp1(t_d, x_d(3,:)-x_d(4,:), t);
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
    % [t, i_L_2, i_L_2_p]
    filename = ['fwr2_solution_' 'i_L_2' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    t        = data(2:end,1);
    i_L_2    = data(2:end,2);
    i_L_2_p  = data(2:end,3);

    filename = ['fwr2_solution_' 'v_34' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    v_34     = data(2:end,2);
    v_34_p   = data(2:end,3);

    filename = ['fwr2_solution_' 'phi_1' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    phi_1    = data(2:end,2);
    phi_1_p  = data(2:end,3);

    filename = ['fwr2_solution_' 'phi_2' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    phi_2    = data(2:end,2);
    phi_2_p  = data(2:end,3);

    filename = ['fwr2_solution_' 'phi_3' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    phi_3    = data(2:end,2);
    phi_3_p  = data(2:end,3);

    filename = ['fwr2_solution_' 'phi_4' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    phi_4    = data(2:end,2);
    phi_4_p  = data(2:end,3);

    filename = ['fwr2_solution_' 'i_L_1' '_tol=' num2str(tol) '.csv'];
    data     = csvread(filename);
    i_L_1    = data(2:end,2);
    i_L_1_p  = data(2:end,3);

    xt_Q   = [v_34'; i_L_2'];
    xt_Q_p = [v_34_p'; i_L_2_p'];
    xt_P   = [i_L_1'];
    xt_P_p = [i_L_1_p'];
    xb_Q   = [phi_1'];
    xb_Q_p = [phi_1_p'];
    xb_P   = [phi_2'; phi_4'];
    xb_P_p = [phi_2_p'; phi_4_p'];
    xb     = [xt_P; xb_Q; xb_P];
    xb_p   = [xt_P_p; xb_Q_p; xb_P_p];
    xt_P_r = NaN(size(xt_P_p));
    xb_Q_r = NaN(size(xb_Q_p));
    xb_P_r = NaN(size(xb_P_p));
    xb_r   = [xt_P_r; xb_Q_r; xb_P_r];
    e_s    = NaN(1,length(t));
    eh     = NaN(1,length(t));
    eb     = NaN(1,length(t));
    for it = 1:length(t)
        if (mod(it, 1000) == 0)
            fprintf('\ntimestep no.: %i/%i\n', it, length(t));
        end

        % start with initial condition or previous timestep
        if it == 1
            xb_r(:,it) = fsolve(@(x) g_i(xt_Q_p(:,it), x(1), x(2), x(3:4), t(it)), [xt_P_0; xb_Q_0; xb_P_0], options);
        else
            xb_r(:,it) = fsolve(@(x) g_i(xt_Q_p(:,it), x(1), x(2), x(3:4), t(it)), xb_r(:,it-1), options);
        end

        e_s(it) = norm(g_i(xt_Q(:,it), xt_P(:,it), xb_Q(:,it), xb_P(:,it), t(it)));
        eh(it)  = norm(g_i(xt_Q_p(:,it), xt_P_p(:,it), xb_Q_p(:,it), xb_P_p(:,it), t(it)));
        eb(it)  = norm(g_i(xt_Q_p(:,it), xb_r(1,it), xb_r(2,it), xb_r(3:4,it), t(it)));
    end

    max(e_s)
    max(eh)
    max(eb)

    semilogy(t, e_s, t, eh, t, eb);

    return

    filename = 'fwr2_solution.csv';
    csvwrite(filename, [t'; xt_Q; xt_Q_p; xb; xb_p; xb_r]');

    filename = 'fwr2_consistency.csv';
    csvwrite(filename, [t'; e_s; eh; eb]');
end
