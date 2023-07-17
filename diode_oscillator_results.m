problem  = 2; % 1, 2

plot_input       = 0; % 0, 1
plot_convergence = 0; % 0, 1
plot_solution    = 0; % 0, 1
plot_consistency = 1; % 0, 1

% only relevant for input/convergence
variable = 'phi_3'; % phi_3, i_L, phi_2

% only relevant for solution
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
        filename = ['do1_doe_i_t_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t      = dlmread(filename);

        filename = ['do1_doe_param_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param    = dlmread(filename);

        filename         = ['do1_doe_lognoisevariance_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance = dlmread(filename);

        filename = ['do1_doe_N_t=' num2str(N_t) '.dat'];
    elseif (problem == 2)
        filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        i_t      = dlmread(filename);

        filename = ['do2_doe_param_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        param    = dlmread(filename);

        filename         = ['do2_doe_lognoisevariance_' variable '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
        lognoisevariance = dlmread(filename);

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
    elseif (problem == 2)
        phi_2 = data(:,5);
        phi_1 = phi_2 + R*i_s(t);
    end

    i_p = find((L == L_p) & (C == C_p));
    X_p = [t(i_p) L(i_p) C(i_p)];

    f = stk_model('stk_gausscov_aniso', 3);

    f.lognoisevariance = lognoisevariance;
    f.param            = param;

    phi_3_p = stk_predict(f, X(i_t,:), phi_3(i_t), X_p);
    phi_2_p = stk_predict(f, X(i_t,:), phi_2(i_t), X_p);

    if (problem == 1)
        i_L_p = stk_predict(f, X(i_t,:), i_L(i_t), X_p);

        R       = 500;
        v_s     = @(t) sin(2*pi*300*t);
        phi_2_r = v_s(t(i_p)) - R*i_L_p.mean;
    elseif (problem == 2)
        phi_1_p = stk_predict(f, X(i_t,:), phi_1(i_t), X_p);

        i_sp    = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);
        phi_2_r = phi_3_p.mean + L_p*i_sp(t(i_p));
    end

    if (problem == 1)
        write_dat ('do1_solution_phi_3.dat', t(i_p)', [phi_3(i_p)'; phi_3_p.mean'; (phi_3_p.mean - phi_3(i_p))']);

        write_dat ('do1_solution_phi_2.dat', t(i_p)', [phi_2(i_p)'; phi_2_p.mean'; (phi_2_p.mean - phi_2(i_p))'; ...
                                                       phi_2_r'; (phi_2_r - phi_2(i_p))'; (phi_2_p.mean - phi_2_r)']);

        write_dat ('do1_solution_i_L.dat', t(i_p)', [i_L(i_p)'; i_L_p.mean'; (i_L_p.mean - i_L(i_p))']);
    elseif (problem == 2)
        write_dat ('do2_solution_phi_3.dat', t(i_p)', [phi_3(i_p)'; phi_3_p.mean'; (phi_3_p.mean - phi_3(i_p))']);

        write_dat ('do2_solution_phi_2.dat', t(i_p)', [phi_2(i_p)'; phi_2_p.mean'; (phi_2_p.mean - phi_2(i_p))'; ...
                                                       phi_2_r'; (phi_2_r - phi_2(i_p))'; (phi_2_p.mean - phi_2_r)']);

        write_dat ('do2_solution_phi_1.dat', t(i_p)', [phi_1(i_p)'; phi_1_p.mean']);
    end
end

if (plot_consistency == 1)
    N_t = 1000;

    filename = ['do2_doe_i_t_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
    i_t      = dlmread(filename);

    filename = ['do2_doe_param_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
    param    = dlmread(filename);

    filename         = ['do2_doe_lognoisevariance_' 'i_L' '_tol=' num2str(tol) '_N_t=' num2str(N_t) '.dat'];
    lognoisevariance = dlmread(filename);

    filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];

    data = dlmread(filename);
    data = data(2:end,:);

    t = data(:,1);

    R    = 500;
    L    = 1.7e-3;
    i_s  = @(t) 1e-4*sin(2*pi*200*t);
    i_sp = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);
    i_L  = i_s(t);

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
    % [t phi_2 phi_2_p]
    data    = dlmread('do2_solution_phi_2.dat');
    phi_2_p = data(2:end,3);
    % [t phi_3 phi_3_p]
    data    = dlmread('do2_solution_phi_3.dat');
    phi_3_p = data(2:end,3);
    % [t i_L i_L_p]
    data  = dlmread('do2_solution_i_L.dat');
    i_L_p = data(2:end,3);

    e_1 = phi_1_p - (phi_2_p + R*i_s(t));
    e_2 = phi_2_p - (phi_3_p + L*i_sp(t));
    e_3 = i_L_p - i_s(t);

    e_c = [e_1'; e_2'; e_3'];
    e_c = norm(e_c, 2, 'cols');

    % 8.69e-3
    max(e_c)
end
