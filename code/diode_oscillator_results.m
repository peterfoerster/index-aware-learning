problem = 2; % 1, 2

plot_input       = 0; % 0, 1
plot_convergence = 0; % 0, 1
plot_solution    = 0; % 0, 1
plot_consistency = 0; % 0, 1
plot_euler       = 1; % 0, 1

% only relevant for input/convergence
variable = 'phi_3'; % phi_3, i_L, phi_2

% only relevant for solution/consistency
tol = 1e-3; % 1e-2, 5e-3, 1e-3

if (plot_input == 1)
    TOL = [1e-2, 5e-3, 1e-3];
    N_t = 1000;

    for tol=TOL
        if (problem == 1)
            filename = ['do1_doe_i_t_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
            i_t      = dlmread(filename);

            filename = ['do1_doe_N_t=' num2str(N_t) '.dat'];
        elseif (problem == 2)
            filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
            i_t      = dlmread(filename);

            filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
        end

        data = dlmread(filename);
        data = data(2:end,:);

        t = data(:,1);
        L = data(:,2);
        C = data(:,3);
        X = [t L C];

        if (problem == 1)
            filename = ['do1_input_' variable '_tol=' num2str(tol) '.dat'];
        elseif (problem == 2)
            filename = ['do2_input_' variable '_tol=' num2str(tol) '.dat'];
        end

        fid = fopen(filename, 'w');
        fprintf(fid, 'x  y\n');
        dlmwrite(fid, X(i_t,2:3), 'delimiter', '  ', 'append', 'on');
        fclose(fid);
    end
end

if (plot_convergence == 1)
    tol = [1e-2, 5e-3, 1e-3];
    N_t = 1000;

    N_s = NaN(length(tol), 1);
    for it=1:length(tol)
       if (problem == 1)
            filename = ['do1_doe_i_t_' variable '_tol=' num2str(tol(it)) '_N_t=' num2str(N_t) '.dat'];
            i_t      = dlmread(filename);
        elseif (problem == 2)
            filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol(it)) '_N_t=' num2str(N_t) '.dat'];
            i_t      = dlmread(filename);
        end

        N_s(it) = length(i_t);
    end

    if (problem == 1)
        filename = ['do1_plot_input_' variable '.dat'];
    elseif (problem == 2)
        filename = ['do2_plot_input_' variable '.dat'];
    end

    fid = fopen(filename, 'w');
    fprintf(fid, 'x  y\n');
    dlmwrite(fid, [N_s tol'], 'delimiter', '  ', 'append', 'on');
    fclose(fid);
end

if (plot_solution == 1)
    N_t = 1000;

    % prediction point
    L_p = 1.7e-3;
    C_p = 220e-9;

    if (problem == 1)
        filename  = ['do1_doe_i_t_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_phi_3 = dlmread(filename);

        filename    = ['do1_doe_param_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_phi_3 = dlmread(filename);

        filename               = ['do1_doe_lognoisevariance_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_phi_3 = dlmread(filename);

        filename  = ['do1_doe_i_t_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_phi_2 = dlmread(filename);

        filename    = ['do1_doe_param_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_phi_2 = dlmread(filename);

        filename               = ['do1_doe_lognoisevariance_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_phi_2 = dlmread(filename);

        filename = ['do1_doe_i_t_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_i_L  = dlmread(filename);

        filename  = ['do1_doe_param_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_i_L = dlmread(filename);

        filename             = ['do1_doe_lognoisevariance_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_i_L = dlmread(filename);

        filename = ['do1_doe_i_t_' 'i_V' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_i_V  = dlmread(filename);

        filename  = ['do1_doe_param_' 'i_V' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_i_V = dlmread(filename);

        filename             = ['do1_doe_lognoisevariance_' 'i_V' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_i_V = dlmread(filename);

        filename = ['do1_doe_N_t=' num2str(N_t) '.dat'];
    elseif (problem == 2)
        filename  = ['do2_doe_i_t_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_phi_3 = dlmread(filename);

        filename    = ['do2_doe_param_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_phi_3 = dlmread(filename);

        filename               = ['do2_doe_lognoisevariance_' 'phi_3' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_phi_3 = dlmread(filename);

        filename  = ['do2_doe_i_t_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_phi_2 = dlmread(filename);

        filename    = ['do2_doe_param_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_phi_2 = dlmread(filename);

        filename               = ['do2_doe_lognoisevariance_' 'phi_2' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_phi_2 = dlmread(filename);

        filename  = ['do2_doe_i_t_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t_phi_1 = dlmread(filename);

        filename    = ['do2_doe_param_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param_phi_1 = dlmread(filename);

        filename               = ['do2_doe_lognoisevariance_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance_phi_1 = dlmread(filename);

        filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
    end

    data = dlmread(filename);
    data = data(2:end,:);

    t = data(:,1);
    L = data(:,2);
    C = data(:,3);
    X = [t L C];

    phi_3 = data(:,4);

    R   = 500;
    i_s = @(t) 1e-4*sin(2*pi*200*t);
    if (problem == 1)
        i_L   = data(:,5);
        phi_2 = data(:,6);
        i_V   = -i_L;
    elseif (problem == 2)
        phi_2 = data(:,5);
        phi_1 = phi_2 + R*i_s(t);
    end

    i_p = find((L == L_p) & (C == C_p));
    X_p = [t(i_p) L(i_p) C(i_p)];

    f = stk_model('stk_gausscov_aniso', 3);

    f.lognoisevariance = lognoisevariance_phi_3;
    f.param            = param_phi_3;

    phi_3_p = stk_predict(f, X(i_t_phi_3,:), phi_3(i_t_phi_3), X_p);

    f.lognoisevariance = lognoisevariance_phi_2;
    f.param            = param_phi_2;

    phi_2_p = stk_predict(f, X(i_t_phi_2,:), phi_2(i_t_phi_2), X_p);

    if (problem == 1)
        f.lognoisevariance = lognoisevariance_i_L;
        f.param            = param_i_L;

        i_L_p = stk_predict(f, X(i_t_i_L,:), i_L(i_t_i_L), X_p);

        f.lognoisevariance = lognoisevariance_i_V;
        f.param            = param_i_V;

        i_V_p = stk_predict(f, X(i_t_i_V,:), i_V(i_t_i_V), X_p);

        R       = 500;
        v_s     = @(t) sin(2*pi*300*t);
        phi_2_r = v_s(t(i_p)) - R*i_L_p.mean;
    elseif (problem == 2)
        f.lognoisevariance = lognoisevariance_phi_1;
        f.param            = param_phi_1;

        phi_1_p = stk_predict(f, X(i_t_phi_1,:), phi_1(i_t_phi_1), X_p);

        i_sp    = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);
        phi_2_r = phi_3_p.mean + L_p*i_sp(t(i_p));
    end

    if (problem == 1)
        write_dat ('do1_solution_phi_3.dat', t(i_p)', [phi_3(i_p)'; phi_3_p.mean'; (phi_3_p.mean - phi_3(i_p))']);

        write_dat ('do1_solution_phi_2.dat', t(i_p)', [phi_2(i_p)'; phi_2_p.mean'; (phi_2_p.mean - phi_2(i_p))'; ...
                                                       phi_2_r'; (phi_2_r - phi_2(i_p))'; (phi_2_p.mean - phi_2_r)']);

        write_dat ('do1_solution_i_L.dat', t(i_p)', [i_L(i_p)'; i_L_p.mean'; (i_L_p.mean - i_L(i_p))']);

        write_dat ('do1_solution_i_V.dat', t(i_p)', [i_V(i_p)'; i_V_p.mean'; (i_V_p.mean - i_V(i_p))']);
    elseif (problem == 2)
        write_dat ('do2_solution_phi_3.dat', t(i_p)', [phi_3(i_p)'; phi_3_p.mean'; (phi_3_p.mean - phi_3(i_p))']);

        write_dat ('do2_solution_phi_2.dat', t(i_p)', [phi_2(i_p)'; phi_2_p.mean'; (phi_2_p.mean - phi_2(i_p))'; ...
                                                       phi_2_r'; (phi_2_r - phi_2(i_p))'; (phi_2_p.mean - phi_2_r)']);

        write_dat ('do2_solution_phi_1.dat', t(i_p)', [phi_1(i_p)'; phi_1_p.mean']);
    end
end

if (plot_consistency == 1)
    N_t = 1000;

    R    = 500;
    L    = 1.7e-3;
    v_s  = @(t) sin(2*pi*300*t);
    i_s  = @(t) 1e-4*sin(2*pi*200*t);
    i_sp = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);

    if (problem == 1)
        filename = ['do1_doe_N_t=' num2str(N_t) '.dat'];
        data     = dlmread(filename);
        data     = data(2:end,:);

        filename = ['do1_doe_i_t_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t      = dlmread(filename);

        filename = ['do1_doe_param_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param    = dlmread(filename);

        filename         = ['do1_doe_lognoisevariance_' 'phi_1' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance = dlmread(filename);

        t     = data(:,1);
        phi_1 = v_s(t);

        i_p = 1:N_t;
        t_p = t(i_p);

        f = stk_model('stk_gausscov_aniso', 1);

        f.lognoisevariance = lognoisevariance;
        f.param            = param;

        phi_1_p = stk_predict(f, t(i_t), phi_1(i_t), t_p);

        filename = 'do1_solution_phi_1.dat';
        if (~exist(filename))
            write_dat (filename, t_p', [phi_1(i_p)'; phi_1_p.mean']);
        end

        % [t phi_1 phi_1_p]
        data    = dlmread('do1_solution_phi_1.dat');
        t       = data(2:end,1);
        phi_1   = data(2:end,2);
        phi_1_p = data(2:end,3);
        phi_1_r = v_s(t);
        % [t phi_2 phi_2_p ... phi_2_r]
        data    = dlmread('do1_solution_phi_2.dat');
        phi_2   = data(2:end,2);
        phi_2_p = data(2:end,3);
        phi_2_r = data(2:end,5);
        % [t phi_3 phi_3_p]
        data    = dlmread('do1_solution_phi_3.dat');
        phi_3   = data(2:end,2);
        phi_3_p = data(2:end,3);
        % [t i_L i_L_p]
        data  = dlmread('do1_solution_i_L.dat');
        i_L   = data(2:end,2);
        i_L_p = data(2:end,3);
        % [t i_V i_V_p]
        data  = dlmread('do1_solution_i_V.dat');
        i_V   = data(2:end,2);
        i_V_p = data(2:end,3);
        i_V_r = -i_L_p;

        e_1 = phi_1_p - phi_1;
        e_2 = phi_2_p - (phi_1 - R*i_L_p);
        e_3 = i_V_p + i_L_p;

        e_c = [e_1'; e_2'; e_3'];
        e_c = norm(e_c, 2, 'cols');

        % 4.3e-3
        max(e_c)

        e_1 = phi_1_r - phi_1;
        e_2 = phi_2_r - (phi_1 - R*i_L_p);
        e_3 = i_V_r + i_L_p;

        e_r = [e_1'; e_2'; e_3'];
        e_r = norm(e_r, 2, 'cols');

        % 1.8e-16
        max(e_r)

        filename = 'do1_consistency.dat';
        write_dat (filename, t', [e_c; e_r+1e-18]);
    elseif (problem == 2)
        filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
        data     = dlmread(filename);
        data     = data(2:end,:);

        filename = ['do2_doe_i_t_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t      = dlmread(filename);

        filename = ['do2_doe_param_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param    = dlmread(filename);

        filename         = ['do2_doe_lognoisevariance_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance = dlmread(filename);

        t   = data(:,1);
        i_L = i_s(t);

        i_p = 1:N_t;
        t_p = t(i_p);

        f = stk_model('stk_gausscov_aniso', 1);

        f.lognoisevariance = lognoisevariance;
        f.param            = param;

        i_L_p = stk_predict(f, t(i_t), i_L(i_t), t_p);

        filename = 'do2_solution_i_L.dat';
        if (~exist(filename))
            write_dat (filename, t_p', [i_L(i_p)'; i_L_p.mean']);
        end

        % [t phi_1 phi_1_p]
        data    = dlmread('do2_solution_phi_1.dat');
        t       = data(2:end,1);
        phi_1_p = data(2:end,3);
        % [t phi_2 phi_2_p ... phi_2_r]
        data    = dlmread('do2_solution_phi_2.dat');
        phi_2_p = data(2:end,3);
        phi_2_r = data(2:end,5);
        phi_1_r = phi_2_r + R*i_s(t);
        % [t phi_3 phi_3_p]
        data    = dlmread('do2_solution_phi_3.dat');
        phi_3_p = data(2:end,3);
        % [t i_L i_L_p]
        data  = dlmread('do2_solution_i_L.dat');
        i_L_p = data(2:end,3);
        i_L_r = i_s(t);

        e_1 = phi_1_p - (phi_2_p + R*i_s(t));
        e_2 = phi_2_p - (phi_3_p + L*i_sp(t));
        e_3 = i_L_p - i_s(t);

        e_c = [e_1'; e_2'; e_3'];
        e_c = norm(e_c, 2, 'cols');

        % 1.6e-3
        max(e_c)

        e_1 = phi_1_r - (phi_2_r + R*i_s(t));
        e_2 = phi_2_r - (phi_3_p + L*i_sp(t));
        e_3 = i_L_r - i_s(t);

        e_r = [e_1'; e_2'; e_3'];
        e_r = norm(e_r, 2, 'cols');

        % 1.1e-16
        max(e_r)

        filename = 'do2_consistency.dat';
        write_dat (filename, t', [e_c; e_r+1e-18]);
    end
end

if (plot_euler == 1)
    % solution: 1e-6
    Deltat = 1e-9;

    R = 500;
    G = 1/R;
    L = 1.7e-3;
    C = 220e-9;

    g_D  = @(v_D) 1e-14*(exp(v_D/26e-3) - 1);
    i_s  = @(t) 1e-4*sin(2*pi*200*t);
    i_sp = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);

    M = [zeros(2,4); zeros(1,2) C 0; zeros(1,3) L];
    K = @(v_D) [G -G 0 0; -G G 0 1; 0 0 g_D(v_D) -1; 0 -1 1 0];
    f = @(t) [i_s(t); 0; 0; 0];

    % load predictions
    % [t phi_1 phi_1_p]
    data    = dlmread('do2_solution_phi_1.dat');
    t       = data(2:end,1);
    phi_1   = data(2:end,2);
    phi_1_p = data(2:end,3);
    % [t phi_2 phi_2_p]
    data    = dlmread('do2_solution_phi_2.dat');
    phi_2   = data(2:end,2);
    phi_2_p = data(2:end,3);
    phi_2_r = data(2:end,5);
    % [t phi_3 phi_3_p]
    data    = dlmread('do2_solution_phi_3.dat');
    phi_3   = data(2:end,2);
    phi_3_p = data(2:end,3);
    % [t i_L i_L_p]
    data  = dlmread('do2_solution_i_L.dat');
    i_L   = data(2:end,2);
    i_L_p = data(2:end,3);

    % compute initial values
    t_i = [t(1); t(2:end) - 1e-9];
    phi_1_i = interp1(t, phi_1_p, t_i);
    phi_2_i = interp1(t, phi_2_p, t_i);
    phi_3_i = interp1(t, phi_3_p, t_i);
    i_L_i   = interp1(t, i_L_p, t_i);

    phi_1_c = NaN(size(phi_1_p));
    phi_2_c = NaN(size(phi_2_p));
    phi_3_c = NaN(size(phi_3_p));
    i_L_c   = NaN(size(i_L_p));

    % compute corrected predictions using two implicit euler steps
    for it=1:length(t_i)
        % initial values
        x_0 = [phi_1_i(it); phi_2_i(it); phi_3_i(it); i_L_i(it)];

        % v_D = phi_3 = x_3
        [x_c, t_c] = implicit_euler (t_i(it), t_i(it)+2*Deltat, Deltat, M, @(x) K(x(3)), @(t) f(t), x_0);

        phi_1_c(it) = x_c(1,end);
        phi_2_c(it) = x_c(2,end);
        phi_3_c(it) = x_c(3,end);
        i_L_c(it)   = x_c(4,end);
    end

    write_dat ('do2_euler_phi_3.dat', t', [phi_3'; phi_3_p'; phi_3_c'; (phi_3_p - phi_3)'; (phi_3_c - phi_3)'; ...
                                           (phi_3_p - phi_3_c)']);

    write_dat ('do2_euler_phi_2.dat', t', [phi_2'; phi_2_p'; phi_2_c'; (phi_2_p - phi_2)'; (phi_2_c - phi_2)'; ...
                                          (phi_2_p - phi_2_c)'; (phi_2_p - phi_2_r)'; (phi_2_c - phi_2_r)']);

    % compute consistency error
    e_1 = phi_1_c - (phi_2_c + R*i_s(t));
    e_2 = phi_2_c - (phi_3_c + L*i_sp(t));
    e_3 = i_L_c - i_s(t);

    e_c = [e_1'; e_2'; e_3'];
    e_c = norm(e_c, 2, 'cols');

    % 1.3e-7
    max(e_c)

    filename = 'do2_euler_consistency.dat';
    write_dat (filename, t', e_c);
end
