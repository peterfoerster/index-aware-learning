T_0    = 0;
T      = 3;
Deltat = 0.1;

G_1 = G_2 = G_3 = G_4 = G_5 = G_6 = G_7 = 1;
L_1 = L_2 = L_3 = L_4 = L_5 = L_6 = L_7 = 1;
L_12 = L_34 = L_56 = 0.1;
C_1 = C_2 = C_3 = C_4 = C_5 = C_6 = C_7 = C_8 = C_9 = 1;

v_s_1 = @(t) sin(4*t) + cos(2*t);
v_s_2 = @(t) sin(3*t) + cos(t);

n_v   = 16;
n_R   = 7;
n_L   = 7;
n_C   = 9;
n_V   = 2;

A_R       = zeros(n_v,n_R);
A_R(1,1)  = 1;
A_R(2,1)  = -1;
A_R(8,2)  = 1;
A_R(9,2)  = -1;
A_R(3,3)  = 1;
A_R(4,3)  = -1;
A_R(10,4) = 1;
A_R(11,4) = -1;
A_R(5,5)  = 1;
A_R(6,5)  = -1;
A_R(12,6) = 1;
A_R(13,6) = -1;
A_R(12,7) = 1;
A_R(15,7) = -1;

A_L       = zeros(n_v,n_L);
A_L(2,1)  = 1;
A_L(3,1)  = -1;
A_L(9,2)  = 1;
A_L(10,2) = -1;
A_L(4,3)  = 1;
A_L(5,3)  = -1;
A_L(11,4) = 1;
A_L(12,4) = -1;
A_L(6,5)  = 1;
A_L(7,5)  = -1;
A_L(13,6) = 1;
A_L(14,6) = -1;
A_L(15,7) = 1;
A_L(16,7) = -1;

A_C       = zeros(n_v,n_C);
A_C(3,1)  = 1;
A_C(10,2) = 1;
A_C(3,2)  = -1;
A_C(5,3)  = 1;
A_C(12,4) = 1;
A_C(5,4)  = -1;
A_C(7,5)  = 1;
A_C(14,6) = 1;
A_C(7,6)  = -1;
A_C(10,7) = -1;
A_C(12,8) = -1;
A_C(16,9) = -1;

A_V      = zeros(n_v,n_V);
A_V(1,1) = 1;
A_V(8,2) = 1;

G      = zeros(n_R,n_R);
G(1,1) = G_1;
G(2,2) = G_2;
G(3,3) = G_3;
G(4,4) = G_4;
G(5,5) = G_5;
G(6,6) = G_6;
G(7,7) = G_7;

L      = zeros(n_L,n_L);
L(1,1) = L_1;
L(1,2) = L_12;
L(2,1) = L_12;
L(2,2) = L_2;
L(3,3) = L_3;
L(3,4) = L_34;
L(4,3) = L_34;
L(4,4) = L_4;
L(5,5) = L_5;
L(5,6) = L_56;
L(6,5) = L_56;
L(6,6) = L_6;
L(7,7) = L_7;

C      = zeros(n_C,n_C);
C(1,1) = C_1;
C(2,2) = C_2;
C(3,3) = C_3;
C(4,4) = C_4;
C(5,5) = C_5;
C(6,6) = C_6;
C(7,7) = C_7;
C(8,8) = C_8;
C(9,9) = C_9;

v_s = @(t) [v_s_1(t); v_s_2(t)];

M = [[A_C*C*A_C.', zeros(n_v,n_L), zeros(n_v,n_V)]; [zeros(n_L,n_v), L, zeros(n_L,n_V)]; [zeros(n_V,n_v), zeros(n_V,n_L), zeros(n_V,n_V)];];
K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L,n_L), zeros(n_L,n_V)]; [-A_V.', zeros(n_V,n_L), zeros(n_V,n_V)]];
f = @(t) [zeros(n_v,1); zeros(n_L,1); v_s(t)];

Q        = W        = zeros(n_v+n_L+n_V,11);
Q(1,1)   = W(1,1)   = 1;
Q(2,2)   = W(2,2)   = 1;
Q(4,3)   = W(4,3)   = 1;
Q(6,4)   = W(6,4)   = 1;
Q(8,5)   = W(8,5)   = 1;
Q(9,6)   = W(9,6)   = 1;
Q(11,7)  = W(11,7)  = 1;
Q(13,8)  = W(13,8)  = 1;
Q(15,9)  = W(15,9)  = 1;
Q(24,10) = W(24,10) = 1;
Q(25,11) = W(25,11) = 1;

P        = V        = zeros(n_v+n_L+n_V,14);
P(3,1)   = V(3,1)   = 1;
P(5,2)   = V(5,2)   = 1;
P(7,3)   = V(7,3)   = 1;
P(10,4)  = V(10,4)  = 1;
P(12,5)  = V(12,5)  = 1;
P(14,6)  = V(14,6)  = 1;
P(16,7)  = V(16,7)  = 1;
P(17,8)  = V(17,8)  = 1;
P(18,9)  = V(18,9)  = 1;
P(19,10) = V(19,10) = 1;
P(20,11) = V(20,11) = 1;
P(21,12) = V(21,12) = 1;
P(22,13) = V(22,13) = 1;
P(23,14) = V(23,14) = 1;

Mt   = V.' * M * P;
Kt_P = V.' * K * P;
Kt_Q = V.' * K * Q;
ft   = @(t) V.' * f(t);
Kb_P = W.' * K * P;
Kb_Q = W.' * K * Q;
fb   = @(t) W.' * f(t);

xt_0 = zeros(14,1);
xb_0 = -Kb_Q \ (Kb_P*xt_0 + fb(T_0));

return

%%%%%%%%%%%%
    M = [zeros(2,5); zeros(1,2) C zeros(1,2); zeros(1,3) L 0; zeros(1,5)];
    K = @(v_D) [G -G 0 0 1; -G G 0 1 0; 0 0 g_D(v_D) -1 0; 0 -1 1 0 0; -1 0 0 0 0];
    f = @(t) [0; 0; 0; 0; -v_s(t)];

    % phi_3 = 0, i_L = 0 => phi_1 = v_s, phi_2 = v_s, i_v = 0
    x_0 = [v_s(T_0); v_s(T_0); 0; 0; 0];

    filename = 'do1.dat';
    % if (~exist(filename))
        % v_D = phi_3 = x_3
        % Deltat         = 1e-7;
        [x_ref, t_ref] = implicit_euler (T_0, T, Deltat, M, @(x) K(x(3)), @(t) f(t), x_0);
        % write_dat (filename, t, x);

        % Deltat  = [1e-4 1e-5 1e-6];
        % e_phi_3 = NaN(1, 3);
        % e_i_L   = NaN(1, 3);
        % for it=1:length(Deltat)
        %     [x, t]      = implicit_euler (T_0, T, Deltat(it), M, @(x) K(x(3)), @(t) f(t), x_0);
        %     e_phi_3(it) = norm(x_ref(3,:) - interp1(t, x(3,:), t_ref)) / norm(x_ref(3,:));
        %     e_i_L(it)   = norm(x_ref(4,:) - interp1(t, x(4,:), t_ref)) / norm(x_ref(4,:));
        % end
    % end
% end
