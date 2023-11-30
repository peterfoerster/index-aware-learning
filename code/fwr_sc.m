% syms, sym, simplify, jacobian, inv, null, rank, +, -, *, /, ^, .'
pkg load symbolic

syms v_1 v_2 v_3 v_4 i_L_1 i_L_2 i_V
syms L_1 L_2 L_12 g_D_1 g_D_2 g_D_3 g_D_4 R v_s
syms g_D

v   = [v_1; v_2; v_3; v_4];
i_L = [i_L_1; i_L_2];

n_v = length(v);
n_L = length(i_L);

A_R = [[0, 0, 0, 0, 0]; [1, 0, -1, 0, 0]; [-1, -1, 0, 0, 1]; [0, 0, 1, 1, -1]];
A_L = [[1 0]; [0 1]; [0, 0]; [0, 0]];
A_V = [[1]; [0]; [0]; [0]];

G = [[g_D_1, 0, 0, 0, 0]; [0, g_D_2, 0, 0, 0]; [0, 0, g_D_3, 0, 0]; [0, 0, 0, g_D_4, 0]; [0, 0, 0, 0, 1/R]];
% G = [[g_D, 0, 0, 0, 0]; [0, g_D, 0, 0, 0]; [0, 0, g_D, 0, 0]; [0, 0, 0, g_D, 0]; [0, 0, 0, 0, 1/R]];
L = [[L_1 L_12]; [L_12 L_2]];

M = [[zeros(n_v,n_v), zeros(n_v,n_L), zeros(n_v,1)]; [zeros(n_L,n_v), L, zeros(n_L,1)]; [zeros(1,n_v), zeros(1,n_L), zeros(1,1)];];
K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L,n_L), zeros(n_L,1)]; [-A_V.', zeros(1,n_L), zeros(1,1)]];
f = [zeros(n_v,1); zeros(n_L,1); v_s];

% first basis functions
Q = W = null(M);
P = V = sym([[zeros(4,2)]; [1, 0]; [0, 1]; [0, 0]]);

% first stage
Mt   = V.' * M * P;
Kt_P = V.' * K * P;
Kt_Q = V.' * K * Q;
ft   = V.' * f;
Kb_P = W.' * K * P;
Kb_Q = W.' * K * Q;
fb   = W.' * f;

% second basis functions
Qb = Wb = null(Kb_Q);
