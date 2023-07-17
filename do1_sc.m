% syms, sym, simplify, jacobian, inv, null, rank, +, -, *, /, ^, .'
pkg load symbolic

syms phi_1 phi_2 phi_3 i_L i_V G C L g_D v_s

M = [[0, 0, 0, 0, 0]; [0, 0, 0, 0, 0]; [0, 0, C, 0, 0]; [0, 0, 0, L, 0]; [0, 0, 0, 0, 0]];

K = [[G, -G, 0, 0, 1]; [-G, G, 0, 1, 0]; [0, 0, g_D, -1, 0]; [0, -1, 1, 0, 0]; [-1, 0, 0, 0, 0]];

f = [[0]; [0]; [0]; [0]; [-v_s]];

% first basis functions
Q = W = sym([[1, 0, 0]; [0, 1, 0]; [0, 0, 0]; [0, 0, 0]; [0, 0, 1]]);
P = V = sym([[0, 0]; [0, 0]; [1, 0]; [0, 1]; [0, 0]]);

% first stage
Mt   = V.' * M * P;
Kt_P = V.' * K * P;
Kt_Q = V.' * K * Q;
ft   = V.' * f;
Kb_P = W.' * K * P;
Kb_Q = W.' * K * Q;
fb   = W.' * f;
