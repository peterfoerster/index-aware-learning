variables = ['i_L_1'; 'i_L_2'; 'v_1'; 'v_2'; 'v_3'; 'v_4'; 'i_V'];

t_0       = 0;
t_f       = 5e-2;
Deltat    = 1e-5;
N_t       = 1000;
epsilon_T = 1e-6;
i_t_i     = [round(1:N_t/10:N_t) N_t];
epsilon_t = 1e-16;
w_1 = w_2 = 1;

for iv = 1:length(variables)
    variable = variables(iv,:);
    t        = linspace(t_0, t_f, N_t);

    % all sampling points
    % T_1 = T_2 = T_3 = T_4 = [20:10:100] + 273.15;
    T_1 = T_2 = [20:10:100] + 273.15;

    % initial points
    N_t_i = length(i_t_i);
    % t_i   = repmat([t_0*ones(2^0, 1); t_f*ones(2^0, 1)], 2^2, 1);
    % T_1_i = repmat([T_1(1)*ones(2^1, 1); T_1(end)*ones(2^1, 1)], 2^1, 1);
    % T_2_i = repmat([T_2(1)*ones(2^2, 1); T_2(end)*ones(2^2, 1)], 2^0, 1);
    t_i   = repmat(t(i_t_i)', 2^2, 1);
    T_1_i = repmat([T_1(1)*ones(N_t_i*1, 1); T_1(end)*ones(N_t_i*1, 1)], 2^1, 1);
    T_2_i = repmat([T_2(1)*ones(N_t_i*2, 1); T_2(end)*ones(N_t_i*2, 1)], 2^0, 1);
    % T_3_i = repmat([T_3(1)*ones(2^3, 1); T_3(end)*ones(2^3, 1)], 2^1, 1);
    % T_4_i = repmat([T_4(1)*ones(2^4, 1); T_4(end)*ones(2^4, 1)], 2^0, 1);

    % define grid (becomes very large very quickly: true optimization worthwhile)
    % N_t   = round((t_f - t_0) / Deltat) + 1;
    N_T_1 = length(T_1);
    N_T_2 = length(T_2);
    % N_T_3 = length(T_3);
    % N_T_4 = length(T_4);
    X     = NaN(N_t*N_T_1*N_T_2, 3);
    y     = NaN(N_t*N_T_1*N_T_2, 1);
    for i_2 = 1:N_T_2
        for i_1 = 1:N_T_1
            filename = ['fwr_Deltat=' num2str(Deltat) '_T_1=' num2str(T_1(i_1)) '_T_2=' num2str(T_2(i_2)) '_T_3=' num2str(T_2(i_2)) '_T_4=' num2str(T_1(i_1)) '.dat'];
            data     = dlmread(filename);
            data     = data(2:end,:);
            t_d      = data(:,1)';
            x_d      = data(:,2:end)';

            i_t                        = (i_1 - 1) + (i_2 - 1)*N_T_1;
            X(i_t*N_t+1:(i_t+1)*N_t,:) = [t', T_1(i_1)*ones(N_t,1), T_2(i_2)*ones(N_t,1)];

            if (strcmp(variable, 'i_L_1'))
                x = interp1(t_d, x_d(5,:), t);
            elseif (strcmp(variable, 'i_L_2'))
                x = interp1(t_d, x_d(6,:), t);
            else
                variable = variable(1:3);
                if (strcmp(variable, 'v_1'))
                    x = interp1(t_d, x_d(1,:), t);
                elseif (strcmp(variable, 'v_2'))
                    x = interp1(t_d, x_d(2,:), t);
                elseif (strcmp(variable, 'v_3'))
                    x = interp1(t_d, x_d(3,:), t);
                elseif (strcmp(variable, 'v_4'))
                    x = interp1(t_d, x_d(4,:), t);
                elseif (strcmp(variable, 'i_V'))
                    x = interp1(t_d, x_d(7,:), t);
                end
            end
            y(i_t*N_t+1:(i_t+1)*N_t) = x';

            % check if solution is parameter dependent
            if (i_1 == i_2 == 1)
                xh = x;
            else
                delta_T = norm(x - xh) / norm(xh);
                xh = x;
            end
        end
    end

    % only learn time dependence if solution is not parameter dependent
    if (delta_T < epsilon_T)
        fprintf('\n%s does not depend on parameters\n', variable);

        % only consider one set of temperatures
        t_i = t_i(1:N_t_i);
        X   = X(1:N_t,1);
        y   = y(1:N_t);
        t   = X(1:N_t,1);
        I   = NaN(length(t_i), 1);
        for i = 1:length(I)
            i_t  = (t_i(i) >= t - epsilon_t & t_i(i) <= t + epsilon_t);
            I(i) = find(i_t);
        end

        GP = stk_model('stk_gausscov_aniso', 1);
    else
        t   = X(:,1);
        T_1 = X(:,2);
        T_2 = X(:,3);
        % T_3 = X(:,4);
        % T_4 = X(:,5);
        I   = NaN(length(t_i), 1);
        for i = 1:length(I)
            i_t  = (t_i(i) >= t - epsilon_t & t_i(i) <= t + epsilon_t);
            i_1  = (T_1_i(i) == T_1);
            i_2  = (T_2_i(i) == T_2);
            % i_3  = (T_3_i(i) == T_3);
            % i_4  = (T_4_i(i) == T_4);
            I(i) = find(i_t & i_1 & i_2);
        end

        GP = stk_model('stk_gausscov_aniso', 3);
    end

    GP.lognoisevariance = NaN;
    tic;
    [theta, sigma] = stk_param_init(GP, X(I,:), y(I));
    fprintf('\nstk_param_init: %d min\n', toc/60);

    tic;
    [GP.param, GP.lognoisevariance] = stk_param_estim(GP, X(I,:), y(I), theta, sigma);
    fprintf('\nstk_param_estim: %d min\n', toc/60);

    tic;
    y_p = stk_predict(GP, X(I,:), y(I), X);
    fprintf('\nstk_predict: %d min\n', toc/60);

    e_p = 1;
    for tol = [1e-2 5e-3 1e-3] % 1e-2, 5e-3, 1e-3
        fprintf('\ntol: %g\n', tol);
        filename = ['fwr_I_' variable '_tol=' num2str(tol) '.dat'];
        if (exist(filename))
            I = dlmread(filename);

            filename = ['fwr_param_' variable '_tol=' num2str(tol) '.dat'];
            GP.param = dlmread(filename);

            filename            = ['fwr_lognoisevariance_' variable '_tol=' num2str(tol) '.dat'];
            GP.lognoisevariance = dlmread(filename);

            tic;
            [GP.param, GP.lognoisevariance] = stk_param_estim(GP, X(I,:), y(I), GP.param, GP.lognoisevariance);
            fprintf('\nstk_param_estim: %d min\n', toc/60);

            tic;
            y_p = stk_predict(GP, X(I,:), y(I), X);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p = norm(y_p.mean - y) / norm(y);
            fprintf('\ne_p: %g\n', e_p);
            continue;
        end

        while (e_p > tol)
            if (delta_T >= epsilon_T)
                fprintf('\niteration no: %i\n', length(I) - N_t_i*2^2 + 1);
            else
                fprintf('\niteration no: %i\n', length(I) - N_t_i + 1);
            end

            % avoid duplicate samples
            y_p.var(I) = zeros(length(I), 1);

            % only if solution is parameter dependent
            if (delta_T >= epsilon_T)
                % scale already sampled temperatures
                T_1_doe = unique(T_1(I));
                for i = 1:length(T_1_doe)
                    i_1          = find(T_1 == T_1_doe(i));
                    y_p.var(i_1) = w_1 * y_p.var(i_1);
                end

                T_2_doe = unique(T_2(I));
                for i = 1:length(T_2_doe)
                    i_2          = find(T_2 == T_2_doe(i));
                    y_p.var(i_2) = w_2 * y_p.var(i_2);
                end
            end

            % maximize variance (bound optimization, derivatives of variance w.r.t. inputs)
            [V_max, i_max] = max(y_p.var);
            I              = [I; i_max];

            tic;
            [GP.param, GP.lognoisevariance] = stk_param_estim(GP, X(I,:), y(I), GP.param, GP.lognoisevariance);
            fprintf('\nstk_param_estim: %d min\n', toc/60);

            tic;
            y_p = stk_predict(GP, X(I,:), y(I), X);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p = norm(y_p.mean - y) / norm(y);
            fprintf('\ne_p: %g\n', e_p);
        end

        filename = ['fwr_I_' variable '_tol=' num2str(tol) '.dat'];
        dlmwrite(filename, I);

        filename = ['fwr_param_' variable '_tol=' num2str(tol) '.dat'];
        dlmwrite(filename, GP.param);

        filename = ['fwr_lognoisevariance_' variable '_tol=' num2str(tol) '.dat'];
        dlmwrite(filename, GP.lognoisevariance);
    end
end
