problem = '2.2'; % 1.1, 1.2, 2.1, 2.2

T_0    = 0;
T      = 1e-2;
Deltat = 1e-6;

R = 500;
G = 1/R;
L = 2e-3;
C = 100e-9;

g_D  = @(v_D) 1e-14*(exp(v_D/26e-3) - 1);

v_s  = @(t) sin(2*pi*300*t); % f = 300
i_s  = @(t) 1e-4*sin(2*pi*200*t); % f = 200
i_sp = @(t) 1e-4*2*pi*200*cos(2*pi*200*t);

% index 1 problem: M x' + K(x) x = f(t)
if (strcmp(problem, '1.1'))
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
end

% decoupled index 1 problem: M x' + K(x) x = f(t)
if (strcmp(problem, '1.2'))
    % index 1: L and C relevant
    for L = [1e-3:0.1e-3:3e-3]
        for C = [100e-9:10e-9:300e-9]
            M = [C 0; 0 L];
            K = @(v_D) [g_D(v_D) -1; 1 1/G];
            f = @(t) [0; v_s(t)];

            x_0 = [0; 0];

            filename = ['do1_L=' num2str(L) '_C=' num2str(C) '.dat'];
            if (~exist(filename))
                % v_D = phi_3 = x_1
                [x, t] = implicit_euler (T_0, T, Deltat, M, @(x) K(x(1)), @(t) f(t), x_0);

                xt = x;
                % algebraic DOFs (i_L = x_2)
                N_t = length(t);
                xb  = [v_s(t); v_s(t) - xt(2,:)/G; -xt(2,:)];

                write_dat (filename, t, [xt; xb]);
            end
        end
    end
end

% index 2 problem: M x' + K(x) x = f(t)
if (strcmp(problem, '2.1'))
    M = [zeros(2,4); zeros(1,2) C 0; zeros(1,3) L];
    K = @(v_D) [G -G 0 0; -G G 0 1; 0 0 g_D(v_D) -1; 0 -1 1 0];
    f = @(t) [i_s(t); 0; 0; 0];

    % phi_3 = 0 => phi_2 = -L i_s' = 0, i_L = i_s, phi_1 = i_s/G + phi_2 = i_s/G
    x_0 = [i_s(T_0)/G; 0; 0; i_s(T_0)];

    filename = 'do2.dat';
    % if (~exist(filename))
        % v_D = phi_3 = x_3
        % Deltat         = 1e-7;
        [x_ref, t_ref] = implicit_euler (T_0, T, Deltat, M, @(x) K(x(3)), @(t) f(t), x_0);
        % write_dat (filename, t, x);

        % Deltat  = [1e-4 1e-5 1e-6];
        % e_phi_3 = NaN(1, 3);
        % for it=1:length(Deltat)
        %     [x, t]      = implicit_euler (T_0, T, Deltat(it), M, @(x) K(x(3)), @(t) f(t), x_0);
        %     % 5e-2, 5e-3, 5e-4
        %     e_phi_3(it) = norm(x_ref(3,:) - interp1(t, x(3,:), t_ref)) / norm(x_ref(3,:));
        % end
    % end
end

% decoupled index 2 problem: M x' + K(x) x = f(t)
if (strcmp(problem, '2.2'))
    % index 2: only C relevant
    for C = [100e-9:10e-9:300e-9]
        M = C;
        K = @(v_D) g_D(v_D);
        f = @(t) i_s(t);

        x_0 = 0;

        % v_D = phi_3 = x
        [x, t] = implicit_euler (T_0, T, Deltat, M, @(x) K(x), @(t) f(t), x_0);

        xt = x;
        for L = [1e-3:0.1e-3:3e-3]
            filename = ['do2_L=' num2str(L) '_C=' num2str(C) '.dat'];

            % algebraic DOFs (phi_2 = phi_3 + L*i_s', phi_1 = phi_2 + i_s/G, i_L = i_s)
            xb = [xt + L*i_sp(t) + i_s(t)/G; xt + L*i_sp(t); i_s(t)];

            write_dat (filename, t, [xt; xb]);
        end
    end
end
