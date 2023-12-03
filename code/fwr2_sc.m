% syms, sym, simplify, jacobian, inv, null, rank, +, -, *, /, ^, .'
pkg load symbolic

syms phi_1 phi_2 phi_3 phi_4 i_L_1 i_L_2
syms L_1 L_2 L_12 g_D_1 g_D_2 g_D_3 g_D_4 R i_s i_sp

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

M = [[zeros(n_phi, n_phi), zeros(n_phi, n_L), zeros(n_phi, n_V)]; [zeros(n_L, n_phi), L, zeros(n_L, n_V)]; [zeros(n_V, n_phi), zeros(n_V, n_L), zeros(n_V, n_V)]];
% K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]; [-A_V.', zeros(n_V, n_L), zeros(n_V, n_V)]];
K = [[A_R*G*A_R.', A_L, A_V]; [-A_L.', zeros(n_L, n_L), zeros(n_L, n_V)]];
f = [A_I*i_s; zeros(n_L, 1); zeros(n_V, 1)];

% first basis functions
Q = W = null(M);
P = V = sym([[zeros(4, 2)]; [1, 0]; [0, 1]]);

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
Pb = Vb = sym([[0, 0, 0]; [1, 0, 0]; [0, 1, 0]; [0, 0, 1]]);

Qt = null(Wb.' * Kb_P);
Pt = sym([1; 0]);

ft_P  = inv(Wb.' * Kb_P * Pt) * Wb.' * fb;
ft_Pp = -i_sp;

Kb   = -inv(Vb.' * Kb_Q * Pb) * Vb.' * Kb_P * Qt;
fb_P = -inv(Vb.' * Kb_Q * Pb) * (Vb.' * Kb_P * Pt + Vb.' * fb);
fh   = Mt * Pt * ft_Pp + Kb * Pb * ft_P + Kt_Q * Pb * fb_P + ft;
return

% (Mt * Qt).';
Wt = sym([[0]; [1]])
Vt = sym([[1]; [0]])

% second stage
Mt_2  = Vt.' * Mt * Qt
Kt_P2 = Vt.' * (Mt * Qtp + Kt_P * Qt)
Kt_Q2 = Vt.' * Kt_Q * Qb
ft_2  = Vt.' * (ft - Mt * Pt * xt_Pp - (Mt * Ptp + Kt_P * Pt) * xt_P)
Kb_P2 = Wt.' * (Mt * Qtp + Kt_P * Qt)
Kb_Q2 = Wt.' * Kt_Q * Qb
fb_2  = Wt.' * (ft - Mt * Pt * xt_Pp - (Mt * Ptp + Kt_P * Pt) * xt_P)
