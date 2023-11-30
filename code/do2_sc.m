% syms, sym, simplify, jacobian, inv, null, rank, +, -, *, /, ^, .'
pkg load symbolic

syms phi_1 phi_2 phi_3 i_L G C L g_D i_s i_sp

M = [[0, 0, 0, 0]; [0, 0, 0, 0]; [0, 0, C, 0]; [0, 0, 0, L]];
K = [[G, -G, 0, 0]; [-G, G, 0, 1]; [0, 0, g_D, -1]; [0, -1, 1, 0]];
f = [[i_s]; [0]; [0]; [0]];

% first basis functions
Q = W = sym([[1, 0]; [0, 1]; [0, 0]; [0, 0]]);
P = V = sym([[0, 0]; [0, 0]; [1, 0]; [0, 1]]);

% first stage
Mt   = V.' * M * P;
Kt_P = V.' * K * P;
Kt_Q = V.' * K * Q;
ft   = V.' * f;
Kb_P = W.' * K * P;
Kb_Q = W.' * K * Q;
fb   = W.' * f;

% second basis functions
Qb = Wb = sym([[1]; [1]])

% Wb.' * Kb_P
Qt  = sym([[1]; [0]])
Qtp = sym([[0]; [0]]);
Pt  = sym([[0]; [1]])
Ptp = sym([[0]; [0]]);

% xt_P = inv(Wb.' * Kb_P * Pt) * Wb.' * fb;
xt_P  = i_s;
xt_Pp = i_sp;

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
