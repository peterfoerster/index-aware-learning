problem  = 1; % 1, 2
variable = 'i_V'; % phi_3, i_L, phi_2, phi_1, i_V
tol      = [1e-2 5e-3 1e-3]; % 1e-2, 5e-3, 1e-3

N_t = 1000;
T_0 = 0;
T   = 1e-2;

% all available data
L = [1e-3:0.1e-3:3e-3];
C = [100e-9:10e-9:300e-9];
t = linspace(T_0, T, N_t);

% initial points
t_i = [t(1); t(end); t(1); t(end); t(1); t(end); t(1); t(end)];
L_i = [L(1); L(1); L(end); L(end); L(1); L(1); L(end); L(end)];
C_i = [C(1); C(1); C(1); C(1); C(end); C(end); C(end); C(end)];

if (problem == 1)
    filename = ['do1_doe_N_t=' num2str(N_t) '.dat'];
elseif (problem == 2)
    filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
end
data = dlmread(filename);
data = data(2:end,:);

t_d = data(:,1);
L_d = data(:,2);
C_d = data(:,3);
X   = [t_d L_d C_d];

phi_3_d = data(:,4);

R   = 500;
i_s = @(t) 1e-4*sin(2*pi*200*t);
if (problem == 1)
    i_L_d   = data(:,5);
    phi_2_d = data(:,6);
    i_V_d   = -i_L_d;
elseif (problem == 2)
    phi_2_d = data(:,5);
    phi_1_d = phi_2_d + R*i_s(t_d);
end

if (strcmp(variable, 'phi_3'))
    y = phi_3_d;
elseif (strcmp(variable, 'i_L'))
    y = i_L_d;
elseif (strcmp(variable, 'phi_2'))
    y = phi_2_d;
elseif (strcmp(variable, 'phi_1'))
    y = phi_1_d;
elseif (strcmp(variable, 'i_V'))
    y = i_V_d;
end

i_i = NaN(length(t_i), 1);
for i=1:length(i_i)
    i_i(i) = find((t_i(i) == t_d) & (L_i(i) == L_d) & (C_i(i) == C_d));
end

f = stk_model('stk_gausscov_aniso', 3);

f.lognoisevariance = NaN;
% f.lognoisevariance = -Inf;
tic;
[theta, sigma] = stk_param_init(f, X(i_i,:), y(i_i));
fprintf('\nstk_param_init: %d min\n', toc/60);

tic;
[f.param, f.lognoisevariance] = stk_param_estim(f, X(i_i,:), y(i_i), theta, sigma);
% [f.param] = stk_param_estim(f, X, y, theta);
fprintf('\nstk_param_estim: %d min\n', toc/60);

tic;
y_p = stk_predict(f, X(i_i,:), y(i_i), X);
fprintf('\nstk_predict: %d min\n', toc/60);

i_t = i_i;
e_p = 1;
for i_tol=1:length(tol)
    fprintf('\ntol: %g\n', tol(i_tol));
    if (problem == 1)
        filename = ['do1_doe_i_t_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        if (exist(filename))
            i_t = dlmread(filename);

            filename = ['do1_doe_param_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
            f.param = dlmread(filename);

            filename = ['do1_doe_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
            f.lognoisevariance = dlmread(filename);

            tic;
            [f.param, f.lognoisevariance] = stk_param_estim(f, X(i_t,:), y(i_t), f.param, f.lognoisevariance);
            fprintf('\nstk_param_estim: %d min\n', toc/60);

            tic;
            y_p = stk_predict(f, X(i_t,:), y(i_t), X);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p = norm(y_p.mean - y) / norm(y);

            fprintf('\ne_p: %g\n', e_p);
        end
    elseif (problem == 2)
        filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        if (exist(filename))
            i_t = dlmread(filename);

            filename = ['do2_doe_param_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
            f.param = dlmread(filename);

            filename = ['do2_doe_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
            f.lognoisevariance = dlmread(filename);

            tic;
            [f.param, f.lognoisevariance] = stk_param_estim(f, X(i_t,:), y(i_t), f.param, f.lognoisevariance);
            fprintf('\nstk_param_estim: %d min\n', toc/60);

            tic;
            y_p = stk_predict(f, X(i_t,:), y(i_t), X);
            fprintf('\nstk_predict: %d min\n', toc/60);

            e_p = norm(y_p.mean - y) / norm(y);

            fprintf('\ne_p: %g\n', e_p);
        end
    end

    while (e_p > tol(i_tol))
        fprintf('\niteration no: %i\n', length(i_t) - length(i_i));

        % avoid duplicate samples
        y_p.var(i_t) = zeros(length(i_t), 1);

        % maximize variance (bound optimization, derivatives of variance w.r.t. inputs)
        [V_max, i_max] = max(y_p.var);
        i_t            = [i_t; i_max];

        tic;
        [f.param, f.lognoisevariance] = stk_param_estim(f, X(i_t,:), y(i_t), f.param, f.lognoisevariance);
        fprintf('\nstk_param_estim: %d min\n', toc/60);

        tic;
        y_p = stk_predict(f, X(i_t,:), y(i_t), X);
        fprintf('\nstk_predict: %d min\n', toc/60);

        e_p = norm(y_p.mean - y) / norm(y);

        fprintf('\ne_p: %g\n', e_p);
    end

    if (problem == 1)
        filename = ['do1_doe_i_t_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, i_t);

        filename = ['do1_doe_param_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.param);

        filename = ['do1_doe_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.lognoisevariance);
    elseif (problem == 2)
        filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, i_t);

        filename = ['do2_doe_param_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.param);

        filename = ['do2_doe_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.lognoisevariance);
    end
end
