N_t = 1000;

T_0 = 0;
T   = 1e-2;

% all available data
L = [1e-3:0.1e-3:3e-3];
C = [100e-9:10e-9:300e-9];
t = [linspace(T_0, T, N_t)]';

x_t   = [];
phi_3 = [];
phi_2 = [];
x_L   = [];
x_C   = [];
for iL=1:length(L)
    for iC=1:length(C)
        filename = ['do2_L=' num2str(L(iL)) '_C=' num2str(C(iC)) '.dat'];
        data     = dlmread(filename);
        data     = data(2:end,:);
        t_data   = data(:,1);
        x_t      = [x_t; t];

        phi_3 = [phi_3; interp1(t_data, data(:,2), t)];
        phi_2 = [phi_2; interp1(t_data, data(:,4), t)];

        x_L = [x_L; L(iL)*ones(N_t, 1)];
        x_C = [x_C; C(iC)*ones(N_t, 1)];
    end
end

filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
write_dat (filename, x_t', [x_L'; x_C'; phi_3'; phi_2']);
