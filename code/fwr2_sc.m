% syms, sym, simplify, jacobian, inv, null, rank, +, -, *, /, ^, .'
pkg load symbolic

syms phi_1 phi_2 phi_3 phi_4 i_L_1 i_L_2
syms C L_1 L_2 L_12 g_D_1 g_D_2 g_D_3 g_D_4 R i_s i_sp

phi = [phi_1; phi_2; phi_3; phi_4];
i_L = [i_L_1; i_L_2];
i_V = [];

n_phi = length(phi);
n_L   = length(i_L);
n_V   = length(i_V);

A_R = [[0, 0, 0, 0, 0]; [1, 0, -1, 0, 0]; [-1, -1, 0, 0, 1]; [0, 0, 1, 1, -1]];
A_L = [[1, 0]; [0, -1]; [0, 0]; [0, 0]];
A_C = [0; 0; 1; -1];
A_V = [];
A_I = [-1; 0; 0; 0];

G = [[g_D_1, 0, 0, 0, 0]; [0, g_D_2, 0, 0, 0]; [0, 0, g_D_3, 0, 0]; [0, 0, 0, g_D_4, 0]; [0, 0, 0, 0, 1/R]];
L = [[L_1, L_12]; [L_12, L_2]];

M = [[A_C*C*A_C.', zeros(n_phi, n_L), zeros(n_phi, n_V)]; [zeros(n_L, n_phi), L, zeros(n_L, n_V)]; [zeros(n_V, n_phi), zeros(n_V, n_L), zeros(n_V, n_V)]];
% K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]; [-A_V.', zeros(n_V, n_L), zeros(n_V, n_V)]];
K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]];
f = [A_I*i_s; zeros(n_L, 1); zeros(n_V, 1)];

% first basis functions
Q = W = null(M);
P = V = sym([[zeros(2, 3)]; [1, 0, 0]; [0, 0, 0]; [0, 1, 0]; [0, 0, 1]]);

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
Pb = Vb = sym([[0, 0]; [1, 0]; [0, 1]]);

Qt = null(Wb.' * Kb_P);
Pt = sym([0; 1; 0]);

Wt = null((Mt * Qt).');
Wt = [0; -L_2; L_12];
Vt = [[1, 0]; [0, L_2]; [0 L_12]];

% alternative ending
Wh = null((Kt_Q * Qb).');
Vh = sym([0; 1; 0]);

% second stage
Kb_Q2 = Wt.' * Kt_Q * Qb
Kb_Q2 = Vh.' * Kt_Q * Qb
