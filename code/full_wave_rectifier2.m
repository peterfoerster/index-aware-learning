samples     = 0; % 0, 1
convergence = 1; % 0, 1

% two Euler steps for ICs
t_0    = -2e-5;
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

% lullo2017: (4) and (7)
L_m  = @(i_L) 3.5e-3 - 1.2e-4*i_L - 1.3e-4*i_L.^2 + 3e-5*i_L.^3 - 1.9e-6*i_L.^4;
L_1  = @(i_L_1) L_m(i_L_1);
L_2  = @(i_L_2) 0.1*L_m(i_L_2);
% coupling coefficient k=0.9
L_12 = @(i_L_1, i_L_2) sqrt(0.9*L_m(i_L_1)*0.1*L_m(i_L_2));

% v_s = @(t) 325*sin(2*pi*50*t);
% v_s = @(t) 12*sin(2*pi*50*t);
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

% first basis functions
Q      = W      = zeros(6,4);
Q(1,1) = W(1,1) = 1;
Q(2,2) = W(2,2) = 1;
Q(3,3) = W(3,3) = 1;
Q(4,4) = W(4,4) = 1;
P      = V      = zeros(6,2);
P(5,1) = V(5,1) = 1;
P(6,2) = V(6,2) = 1;

Mt   = @(x) V.' * M(x)*P;
% Kt_P = @(x, T) V.' * K(x, T)*P;
% Kt_Q = @(x, T) V.' * K(x, T)*Q;
% ft   = @(t) V.' * f(t);
% Kb_P = @(x, T) W.' * K(x, T)*P;
Kt_P = [0, 0; 0, 0];
Kt_Q = [-1, 0, 0, 0; 0, 1, 0, 0];
ft   = [0; 0];
Kb_P = [1, 0; 0, -1; 0, 0; 0, 0];
Kb_Q = @(x, T) W.' * K(x, T)*Q;
fb   = @(t) W.' * f(t);

% second basis functions
Qb      = Wb      = zeros(4,1);
Qb(1,1) = Wb(1,1) = 1;
Pb      = Vb      = zeros(4,3);
Pb(2,1) = Vb(2,1) = 1;
Pb(3,2) = Vb(3,2) = 1;
Pb(4,3) = Vb(4,3) = 1;

Qt = [0; 1];
Pt = [1; 0];

% function handle at x!
Wt = @(x) [-L_2(x(6)); L_12(x(5), x(6))];
Vt = @(x) [L_2(x(6)); L_12(x(5), x(6))];
Wh = [0; 1];
Vh = [1; 0];

T       = [60; 80; 80; 60] + 273.15;
% i_L_2
xt_Q_0  = 0;
xt_P    = @(t) i_s(t);
xt_Pp   = @(t) i_sp(t);
% i_L_1
xt_P_0  = xt_P(t_0);
% phi_1
xb_Q_0  = 0;
% phi_2, phi_3, phi_4
xb_P_0  = [0; 0; 0];

function [g] = g(xt_Q, xt_P, xt_Pp, xb_Q, xb_P, t, T, Deltat, Q, P, Qb, Pb, Wb, Vb, Qt, Pt, Wt, Mt, Kt_P, Kt_Q, ft, Kb_P, Kb_Q, fb, i_s, i_sp)
    xt_P    = -(Wb.' * Kb_P*Pt) \ (Wb.' * fb(t));
    xt_P_po = -(Wb.' * Kb_P*Pt) \ (Wb.' * fb(t+Deltat));

    % approximate derivative by backward difference
    xt_Pp      = (xt_P_po - xt_P) / Deltat;
    Deltaxt_P  = xt_P - i_s(t);
    Deltaxt_Pp = xt_Pp - i_sp(t);

    % x = P xt + Q xb = P (Pt xt_P + Qt xt_Q) + Q (Pb xb_P + Qb xb_Q)
    x      = P*(Pt*xt_P + Qt*xt_Q) + Q*(Pb*xb_P + Qb*xb_Q);

    fb_P   = -(Vb.' * Kb_Q(x, T)*Pb) \ (Vb.' * Kb_P*Pt*xt_P + Vb.' * fb(t));
    fh     = Mt(x)*Pt*xt_Pp + Kt_P*Pt*xt_P + Kt_Q*Pb*fb_P + ft;
    Kbxt_Q = -(Vb.' * Kb_Q(x, T)*Pb) \ (Vb.' * Kb_P*Qt*xt_Q);

    % solve joint nonlinear problem for all algebraic DOFs
    % g_xt_P = Wb.' * Kb_P(x, T)*Pt*xt_P + Wb.' * fb(t);
    g_xb_P = Vb.' * Kb_Q(x, T)*Pb*xb_P + Vb.' * Kb_P*(Qt*xt_Q + Pt*xt_P) + Vb.' * fb(t);
    g_xb_Q = Wt(x).' * Kt_Q*Qb*xb_Q + Wt(x).' * (Kt_P*Qt*xt_Q + Kt_Q*Pb*Kbxt_Q) + Wt(x).' * fh;

    g = [g_xb_P; g_xb_Q];
end

% determine consistent initial conditions
g_2            = @(xt_Q, xt_P, xb_Q, xb_P, t) g(xt_Q, xt_P, xt_Pp, xb_Q, xb_P, t, T, Deltat, Q, P, Qb, Pb, Wb, Vb, Qt, Pt, Wt, Mt, Kt_P, Kt_Q, ft, Kb_P, Kb_Q, fb, i_s, i_sp);
options.TolX   = 1e-16;
options.TolFun = 1e-16;
x_0            = fsolve(@(x) g_2(xt_Q_0, x(1), x(2), x(3:5), t_0), [xt_P_0; xb_Q_0; xb_P_0], options);
xt_P_0         = x_0(1);
xb_Q_0         = x_0(2);
xb_P_0         = x_0(3:5);
x_0            = P*(Pt*xt_P_0 + Qt*xt_Q_0) + Q*(Pb*xb_P_0 + Qb*xb_Q_0);

if (samples)
    % T_1 = T_2 = [20:10:100] + 273.15;
    for T_1 = [20:10:100] + 273.15
        for T_2 = [20:10:100] + 273.15
            filename = ['csv/fwr2/fwr2_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1) '_T_2=' num2str(T_2) '_T_3=' num2str(T_2) '_T_4=' num2str(T_1) '.csv'];
            if (~exist(filename))
                T      = [T_1; T_2; T_2; T_1];
                % M x' + K x + f(t) = 0
                [x, t] = implicit_euler(t_0, t_f, Deltat, @(x) M(x), @(x) K(x, T), @(t) -f(t), x_0);
                csvwrite(filename, [t; x]');
            end
        end
    end
end

if (convergence)
    % Deltath  = 1e-7;
    % Deltath  = 1e-6;
    Deltath  = 1e-5;
    filename = ['fwr2_Deltat=' num2str(Deltath) '.csv'];
    if (~exist(filename))
        % M x' + K x + f(t) = 0
        [xh, th] = implicit_euler(t_0, t_f, Deltath, @(x) M(x), @(x) K(x, T), @(t) -f(t), x_0);
        csvwrite(filename, [th; xh]');
    else
        data = csvread(filename);
        th   = data(:,1)';
        xh   = data(:,2:end)';
    end

    % Deltat  = [1e-3 1e-4 1e-5 1e-6];
    % Deltat  = [1e-3 1e-4 1e-5];
    Deltat  = [1e-3 1e-4];
    e_i_L_2 = NaN(size(Deltat));
    for i_t = 1:length(Deltat)
        filename = ['fwr2_Deltat=' num2str(Deltat(i_t)) '.csv'];
        if (~exist(filename))
            [x, t] = implicit_euler(t_0, t_f, Deltat(i_t), @(x) M(x), @(x) K(x, T), @(t) -f(t), x_0);
            csvwrite(filename, [t; x]');
        else
            data = csvread(filename);
            t    = data(:,1)';
            x    = data(:,2:end)';
        end

        i_L_2        = interp1(t(3:end), x(6,3:end), th(3:end));
        e_i_L_2(i_t) = norm(i_L_2 - xh(6,3:end))/norm(xh(6,3:end));
    end

    e_c = NaN(1,length(th));
    for it = 1:length(th)
        if (mod(it, 1000) == 0)
            fprintf('\ntimestep no.: %i/%i\n', it, length(th));
        end

        e_c(it) = norm(g_2(xh(6,it), xh(5,it), xh(1,it), xh(2:4,it), th(it)));
    end

    max(e_c)
    semilogy(th, e_c, 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('t');
    ylabel('e_c');

    figure;
    plot(th(1:end), xh(6,1:end), 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('t');
    ylabel('i_{L_2}');
end
