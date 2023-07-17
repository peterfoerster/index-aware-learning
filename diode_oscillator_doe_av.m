problem  = 2; % 2
variable = 'i_L'; % i_L
tol      = [1e-3]; % 1e-3

N_t = 1000;
T_0 = 0;
T   = 1e-2;

% all available data
t = linspace(T_0, T, N_t);

% initial points (need more than boundary)
t_i = [t(1); t(500); t(end)];

filename = ['do2_doe_N_t=' num2str(N_t) '.dat'];
data     = dlmread(filename);
data     = data(2:end,:);

t_d = data(1:N_t,1);
X   = [t_d];

R   = 500;
i_s = @(t) 1e-4*sin(2*pi*200*t);
if (strcmp(variable, 'i_L'))
    y = i_s(t_d);
end

i_i = NaN(length(t_i), 1);
for i=1:length(i_i)
    i_i(i) = find((t_i(i) == t_d));
end

f = stk_model('stk_gausscov_aniso', 1);

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
    if (problem == 2)
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

            y_p = stk_predict(f, X(i_t,:), y(i_t), X);
            e_p = norm(y_p.mean - y) / norm(y);
        end
    end

    while (e_p > tol(i_tol))
        fprintf('\niteration no: %i\n', length(i_t) - length(i_i));

        % avoid duplicate samples
        y_p.var(i_t) = zeros(length(i_t), 1);

        % maximize variance
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

    if (problem == 2)
        filename = ['do2_doe_i_t_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, i_t);

        filename = ['do2_doe_param_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.param);

        filename = ['do2_doe_lognoisevariance_' variable '_tol=' num2str(tol(i_tol)) '_N_t=' num2str(N_t) '.dat'];
        dlmwrite(filename, f.lognoisevariance);
    end
end
