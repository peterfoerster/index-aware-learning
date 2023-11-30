samples     = 0; % 0, 1
convergence = 1; % 0, 1

t_0    = 0;
t_f    = 5e-2;
Deltat = 1e-5;
N_t    = 1000;

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

% v_s = @(t) 325*sin(2*pi*50*t);
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
Kt_P = @(x, T) V.' * K(x, T) * P;
Kt_Q = @(x, T) V.' * K(x, T) * Q;
ft   = @(t) V.' * f(t);
Kb_P = @(x, T) W.' * K(x, T) * P;
Kb_Q = @(x, T) W.' * K(x, T) * Q;
fb   = @(t) W.' * f(t);

T    = [25; 25; 25; 25] + 273.15;
xt_0 = zeros(2,1);
xb_0 = zeros(5,1);
g    = @(xt, xb, t) Kb_P(P*xt + Q*xb, T) * xt + Kb_Q(P*xt + Q*xb, T) * xb + fb(t);
x_0  = P*xt_0 + Q*xb_0;

if (samples)
    % T_1 = T_2 = [20:10:100] + 273.15;
    for T_1 = [20:10:100] + 273.15
        for T_2 = [20:10:100] + 273.15
            filename = ['dat/fwr/fwr_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1) '_T_2=' num2str(T_2) '_T_3=' num2str(T_2) '_T_4=' num2str(T_1) '.dat'];
            if (~exist(filename))
                T      = [T_1; T_2; T_2; T_1];
                % [x, t] = trapezoidal_rule_MNA(t_0, t_f, Deltat, M, @(x) K(x, T), @(t) f(t), x_0);
                [x, t] = implicit_euler(t_0, t_f, Deltat, M, @(x) K(x, T), @(t) -f(t), x_0);
                write_dat (filename, t, x);
            end
        end
    end
end

if (convergence)
    % Deltath  = 1e-7;
    % Deltath  = 1e-6;
    Deltath  = 1e-5;
    filename = ['fwr_Deltat=' num2str(Deltath) '.dat'];
    if (~exist(filename))
        % [xh, th] = trapezoidal_rule_MNA(t_0, t_f, Deltath, M, @(x) K(x, T), @(t) f(t), x_0);
        [xh, th] = implicit_euler(t_0, t_f, Deltath, M, @(x) K(x, T), @(t) -f(t), x_0);
        % M x' + K x + f(t) = 0
        % options = odeset('Mass', M);
        % [th, xh] = ode15s(@(t, x) -K(x, T)*x - f(t), linspace(t_0, t_f, N_t), x_0, options);
        % [th, xh] = ode15i(@(t, x, xp) M*xp + K(x, T)*x + f(t), linspace(t_0, t_f, N_t), x_0, x_0);
        write_dat (filename, th, xh);
    else
        data = dlmread(filename);
        data = data(2:end,:);
        th   = data(:,1)';
        xh   = data(:,2:end)';
    end

    % Deltat  = [1e-3 1e-4 1e-5 1e-6];
    Deltat  = [1e-3 1e-4 1e-5];
    e_i_L_1 = NaN(size(Deltat));
    e_i_L_2 = NaN(size(Deltat));
    for i_t = 1:length(Deltat)
        filename = ['fwr_Deltat=' num2str(Deltat(i_t)) '.dat'];
        if (~exist(filename))
            % [x, t] = trapezoidal_rule_MNA(t_0, t_f, Deltat(i_t), M, @(x) K(x, T), @(t) f(t), x_0);
            [x, t] = implicit_euler(t_0, t_f, Deltat(i_t), M, @(x) K(x, T), @(t) -f(t), x_0);
            write_dat (filename, t, x);
        else
            data = dlmread(filename);
            data = data(2:end,:);
            t    = data(:,1)';
            x    = data(:,2:end)';
        end

        i_L_1        = interp1(t, x(5,:), th);
        e_i_L_1(i_t) = norm(i_L_1 - xh(5,:))/norm(xh(5,:));

        i_L_2        = interp1(t, x(6,:), th);
        e_i_L_2(i_t) = norm(i_L_2 - xh(6,:))/norm(xh(6,:));
    end

    e_c = NaN(1,length(th));
    for it = 1:length(th)
        if (mod(it, 1000) == 0)
            fprintf('\ntimestep no.: %i/%i\n', it, length(th));
        end

        e_c(it) = norm(g(xh(5:6,it), [xh(1:4,it); xh(7,it)], th(it)));
    end

    max(e_c)
    semilogy(th, e_c, 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('t');
    ylabel('e_c');

    figure;
    plot(th, xh(5,:), 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('t');
    ylabel('i_{L_1}');

    figure;
    plot(th, xh(6,:), 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('t');
    ylabel('i_{L_2}');
end
