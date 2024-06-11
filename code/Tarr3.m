data = readtable('Copia de cps09mar.xlsx');

% Solo hombres (female == 0)
men_data = data(data.female == 0, :);

%% Pregunta 1

%1
% Selección variables
Y = men_data.union; % variable dependiente
X = [men_data.age, men_data.education, dummyvar(categorical(men_data.race)), dummyvar(categorical(men_data.hisp))];

% Añadir una columna de unos para el término constante
X = [ones(height(men_data),1), X];

% Estima los coeficientes usando MCO
beta_pool = inv(X'*X)*(X'*Y); 

% Calcula los residuos
residuopool = Y - X*beta_pool;

% a) Errores estándar asumiendo homocedasticidad y ausencia de correlación
n = size(X,1); % número de observaciones
k = size(X,2); % número de regresores
V_hom = (residuopool'*residuopool / (n - k)) * inv(X'*X); % matriz de varianzas y covarianzas
err_hom_pool = sqrt(diag(V_hom)); % errores estándar

% b) Errores estándar robustos
D = diag(residuopool.^2);
V_rob = inv(X'*X) * (X'*D*X) * inv(X'*X);
err_rob_pool = sqrt(diag(V_rob));

% Muestra los resultados
disp('Coeficientes del Modelo de Probabilidad Lineal:');
disp(beta_pool);

disp('Errores estándar (homocedasticidad):');
disp(err_hom_pool);

disp('Errores estándar robustos:');
disp(err_rob_pool);

%corregir errores estándar...

%%
% Definir la función de verosimilitud negativa
neg_log_likelihood = @(theta) -sum(Y .* log(normcdf(X * theta)) + (1 - Y) .* log(1 - normcdf(X * theta)));

% Estimación de los coeficientes usando la función de optimización fminsearch
theta_init = zeros(size(X, 2), 1); % Inicialización de los parámetros
options = optimset('Display', 'off'); % Opciones de optimización
[beta_probit, fval, exitflag, output] = fminsearch(neg_log_likelihood, theta_init, options);


% Calcular los errores estándar robustos
residuals = Y - normcdf(X * beta_probit);
D = diag(residuals.^2);
V_rob_probit = inv(X' * X) * (X' * D * X) * inv(X' * X);
std_errors_robust_probit = sqrt(diag(V_rob_probit));

% Mostrar los resultados
disp('Coeficientes del Modelo Probit:');
disp(beta_probit);

disp('Errores estándar robustos del Modelo Probit:');
disp(std_errors_robust_probit);

%Corregir errores estándar...

%% Pregunta 2
%1
% Parámetros dados
rng(123);
beta0 = 1;
beta1 = 2;
beta2 = 5;
alpha0 = -4;
alpha2 = 3;
N = 1000;
alphas = [0.1, 0.5, 1, 5, 10]; % Valores de alpha1 a considerar

% Preparar matrices para almacenar resultados
coeficientes = zeros(length(alphas), 2);
errores_estandar = zeros(length(alphas), 2);

for i = 1:length(alphas)
    alpha1 = alphas(i);
    
    % Generación de datos
    W = normrnd(2, 1, N, 1);  % Genera Wi ~ N(2,1)
    ei = normrnd(0, 1, N, 1); % Genera ei ~ N(0,1)
    ui = normrnd(0, 1, N, 1); % Genera ui ~ N(0,1)
    vi = rand(N, 1);          % Genera νi ~ U[0,1]
    Z = zeros(N, 1);
    Z(vi < 0.8) = 1;          % Zi = 1 si νi < 0.8

    Xi = alpha0 + alpha1 * Z + alpha2 * W + ui; % Xi según ecuación
    Y = beta0 + beta1 * Xi + ei; % Yi según ecuación

    % Estimación usando MCO
    X = [ones(N, 1), Xi]; % Matriz de regresores
    b = (X' * X) \ (X' * Y); % Estimación de coeficientes

    % Cálculo de los errores estándar
    e = Y - X * b; % Residuos
    sigma2 = sum(e.^2) / (N - 2); % Varianza del error
    SE = sqrt(sigma2 * diag(inv(X' * X))); % Errores estándar de los coeficientes

    % Guardar resultados
    coeficientes(i, :) = b';
    errores_estandar(i, :) = SE';
end




% Reporte de resultados
disp('Resultados de la estimación usando MCO para diferentes valores de alpha1:');
for i = 1:length(alphas)
    disp(['alpha1 = ', num2str(alphas(i))]);
    disp(['Coeficientes:']);
    disp(['beta0 = ', num2str(coeficientes(i, 1))]);
    disp(['beta1 = ', num2str(coeficientes(i, 2))]);
    disp(['Errores estándar:']);
    disp(['SE(beta0) = ', num2str(errores_estandar(i, 1))]);
    disp(['SE(beta1) = ', num2str(errores_estandar(i, 2))]);
    disp(' ');
end


%% CREO QUE NO ES LO QUE PIDEN...
% Parámetros dados
beta0 = 1;
beta1 = 2;
beta2 = 5;
alpha0 = -4;
alpha2 = 3;
N = 1000;
alphas = [0.1, 0.5, 1, 5, 10]; % Valores de alpha1 a considerar

% Generación de datos
W = normrnd(2, 1, N, 1);  % Genera Wi ~ N(2,1)
ei = normrnd(0, 1, N, 1); % Genera ei ~ N(0,1)
ui = normrnd(0, 1, N, 1); % Genera ui ~ N(0,1)
vi = rand(N, 1);          % Genera νi ~ U[0,1]
Z = zeros(N, 1);
Z(vi < 0.8) = 1;          % Define Zi = 1 si νi < 0.8

Xi = alpha0 + alphas(1) * Z + alpha2 * W + ui; % Genera Xi según la ecuación dada
Y = beta0 + beta1 * Xi + ei; % Genera Yi según la ecuación dada

% Estimación usando MCO
X = [ones(N, 1), Xi]; % Matriz de regresores
b = (X' * X) \ (X' * Y); % Estimación de coeficientes

% Cálculo de los errores estándar
e = Y - X * b; % Residuos
sigma2 = sum(e.^2) / (N - 2); % Varianza del error
SE = sqrt(sigma2 * diag(inv(X' * X))); % Errores estándar de los coeficientes

% Reporte de resultados
disp('Resultados de la estimación usando MCO:')
disp(['Coeficientes:']);
disp(['beta0 = ', num2str(b(1))]);
disp(['beta1 = ', num2str(b(2))]);
disp([' ']);
disp(['Errores estándar:']);
disp(['SE(beta0) = ', num2str(SE(1))]);
disp(['SE(beta1) = ', num2str(SE(2))]);
%%
% Parámetros dados
beta0_true = 1;
beta1_true = 2;
alpha0 = -4;
alpha2 = 3;
N = 1000;
alphas = [0.1, 0.5, 1, 5, 10]; % Valores de alpha1 a considerar
num_simulations = 1000;

% Preparar matrices para almacenar resultados
beta_hat = zeros(num_simulations, 2); % Almacenar beta0_hat y beta1_hat
bias = zeros(num_simulations, 2); % Almacenar bias de beta0_hat y beta1_hat

for sim = 1:num_simulations
    % Generación de datos
    W = normrnd(2, 1, N, 1);  % Genera Wi ~ N(2,1)
    ei = normrnd(0, 1, N, 1); % Genera ei ~ N(0,1)
    ui = normrnd(0, 1, N, 1); % Genera ui ~ N(0,1)
    vi = rand(N, 1);          % Genera νi ~ U[0,1]
    Z = zeros(N, 1);
    Z(vi < 0.8) = 1;          % Define Zi = 1 si νi < 0.8
    
    % Genera Xi según diferentes valores de alpha1
    for j = 1:length(alphas)
        alpha1 = alphas(j);
        Xi = alpha0 + alpha1 * Z + alpha2 * W + ui; % Genera Xi según la ecuación dada
        Y = beta0_true + beta1_true * Xi + ei; % Genera Yi según la ecuación dada
        
        % Estimación usando MCO
        X = [ones(N, 1), Xi]; % Matriz de regresores
        b = (X' * X) \ (X' * Y); % Estimación de coeficientes
        
        % Guardar estimaciones
        beta_hat(sim, :) = b';
    end
end

% Cálculo del sesgo
bias(:, 1) = beta_hat(:, 1) - beta0_true;
bias(:, 2) = beta_hat(:, 2) - beta1_true;

% Graficar la distribución asintótica de beta_hat
figure;
subplot(1, 2, 1);
histogram(beta_hat(:, 1), 'Normalization', 'pdf');
title('Distribución Asintótica de \beta_0^{MCO}');
xlabel('\beta_0^{MCO}');
ylabel('Densidad');

subplot(1, 2, 2);
histogram(beta_hat(:, 2), 'Normalization', 'pdf');
title('Distribución Asintótica de \beta_1^{MCO}');
xlabel('\beta_1^{MCO}');
ylabel('Densidad');

% Reporte del sesgo promedio
avg_bias = mean(bias);
disp('Sesgo promedio de beta_hat:');
disp(['Bias(beta0_hat) = ', num2str(avg_bias(1))]);
disp(['Bias(beta1_hat) = ', num2str(avg_bias(2))]);

%% 3


% Parámetros dados
beta0_true = 1;
beta1_true = 2;
alpha0 = -4;
alpha2 = 3;
N = 1000;

% Generación de datos
W = normrnd(2, 1, N, 1);  % Genera Wi ~ N(2,1)
ei = normrnd(0, 1, N, 1); % Genera ei ~ N(0,1)
ui = normrnd(0, 1, N, 1); % Genera ui ~ N(0,1)
vi = rand(N, 1);          % Genera νi ~ U[0,1]
Z = zeros(N, 1);
Z(vi < 0.8) = 1;          % Define Zi = 1 si νi < 0.8

% Genera Xi según diferentes valores de alpha1
alpha1 = 1; % Escoge un valor específico para alpha1
Xi = alpha0 + alpha1 * Z + alpha2 * W + ui; % Genera Xi según la ecuación dada
Y = beta0_true + beta1_true * Xi + ei; % Genera Yi según la ecuación dada

% Estimación usando MC2E (IV)
% Primera etapa: Estimar Xi usando Zi como instrumento
X_first_stage = [ones(N, 1), Z];
b_first_stage = (X_first_stage' * X_first_stage) \ (X_first_stage' * Xi);

% Segunda etapa: Estimar Y usando Xi estimado en la primera etapa
X_second_stage = [ones(N, 1), X_first_stage(:, 2)];
b_second_stage = (X_second_stage' * X_second_stage) \ (X_second_stage' * Y);

% Calculo de errores estándar robustos para IV
e_IV = Y - X_second_stage * b_second_stage; % Residuos IV
sigma2_IV = e_IV' * e_IV / N; % Varianza del error IV
SE_IV = sqrt(sigma2_IV * diag(inv(X_second_stage' * X_second_stage))); % Errores estándar IV

% Test F de la primera etapa
RSS_reduced = sum((Xi - X_first_stage * b_first_stage).^2);
RSS_full = sum((Xi - X_first_stage * inv(X_first_stage' * X_first_stage) * X_first_stage' * Xi).^2);
F_statistic = ((RSS_reduced - RSS_full) / 1) / (RSS_full / (N - 2));

% Reporte de resultados
disp('Resultados de la estimación usando MC2E (IV):');
disp(['Coeficientes:']);
disp(['beta0_IV = ', num2str(b_second_stage(1))]);
disp(['beta1_IV = ', num2str(b_second_stage(2))]);
disp([' ']);
disp(['Errores estándar robustos:']);
disp(['SE(beta0_IV) = ', num2str(SE_IV(1))]);
disp(['SE(beta1_IV) = ', num2str(SE_IV(2))]);
disp([' ']);
disp(['Test F de la primera etapa:']);
disp(['Estadístico F = ', num2str(F_statistic)]);
