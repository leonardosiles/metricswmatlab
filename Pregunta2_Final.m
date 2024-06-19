%% Pregunta 2: Tarea 3 (Econometría I)

clear all
clc

rng(123)               % Semilla. 

%% Pregunta 1: W_i no observable.

% Si W_i no es observable, tampoco entrará dentro de nuestro instrumento.
% De esta forma, la primera y segunda etapa serán respectivamente:

for a1 = [0.1, 0.5, 1, 5, 10]

b0 = 1;
b1 = 2;
b2 = 5;

a0 = -4;
a2 = 3; 

N = 1000;

e_i = randn(N,1);
u_i = randn(N,1);
W_i = normrnd(2, 1, N, 1);
v_i = rand(N, 1); 

Z_i = zeros(N,1); 
for i = 1:N
    if v_i(i) < 0.8                       
       Z_i(i) = 1;  
    else             
       Z_i(i) = 0;
    end
end

X_i = a0 + a1*Z_i + a2*W_i + u_i;    % Primera etapa.
Y = b0 + b1*X_i + b2*W_i + e_i;    % Segunda etapa.

X = [ones(N,1), X_i];

beta_gorro = inv(X'*X)*X'*Y;

u_gorro = Y - X*beta_gorro;

K = length(beta_gorro);

% Errores estándares.

s = (u_gorro'*u_gorro)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    % Error estándar homocedástico.
se_r = sqrt(diag(inv(X'*X)*(X'*(diag(u_gorro.^2))*(N/(N-K))*X)*inv(X'*X))); % Error estándar robusto.

Pregunta1 = table(beta_gorro, se, se_r, ...
    'VariableNames', {'Coeficiente', 'SE', 'SE Robusto'});
disp(['Para un alpha_1 = ', num2str(a1), ', los coeficientes y errores estándares son:']);
disp(Pregunta1);
end

%% Pregunta 2: Distribución asintótica por Montecarlo.

% Se repite el mismo procedimiento de la pregunta a), pero considerando las
% B = 1000 simulaciones para beta.

for a1 = [0.1, 0.5, 1, 5, 10]

B = 1000;         % Número de simulaciones. 
bm = NaN(2,B);    % Vector de coeficientes simulados por Montecarlo.

for i = 1:B
    X_i = a0 + a1*Z_i + a2*W_i + u_i;    % Primera etapa.
    Y = b0 + b1*X_i + b2*W_i + e_i;    % Segunda etapa.

    X = [ones(N,1), X_i];   % Primera etapa omitiendo W_i.
    
    m = randi(N,N,1);      % Vector de m pares aleatorios intependientes de tamaño n. 
    bm(:,i) = (X(m,:)'*X(m,:))\(X(m,:)'*Y(m)) ;
end
bm = sort(bm,2);

% Sesgo.
sesgo = mean(bm(2,:))- b1;
disp(['Para un alpha_1 = ', num2str(a1), ', el sesgo de nuestro estimador de MCO será:']);
    disp(sesgo);

% Calcular la densidad estimada
hold on;
[f,xi] = ksdensity(bm(2, :));
plot(xi,f,'LineWidth',2);
hold off;
end

xlabel('Valor');
ylabel('Densidad');
title('Curva de densidad estimada $\hat{\beta}_{1}$', 'Interpreter','latex');
legend('$\alpha_1 = 0.1$', '$\alpha_1 = 0.5$','$\alpha_1 = 1$', ...
    '$\alpha_1 = 5$', '$\alpha_1 = 10$','Interpreter','latex', ...
    'Location', 'northwest');
filename = ['densidad_alpha_p1.eps'];
    print(gcf, filename, '-depsc', '-r300');

%% Pregunta 3: Estimación por MC2E.
clear all   % Se hace clear all. De lo contrario los resultados se traslapan.

rng(123)               % Semilla.

for a1 = [0.1, 0.5, 1, 5, 10]

N = 1000;

b0 = 1;
b1 = 2;
b2 = 5;

a0 = -4;
a2 = 3; 

e_i = randn(N,1);
u_i = randn(N,1);
W_i = normrnd(2, 1, N, 1);
v_i = rand(N, 1); 

Z_i = zeros(N,1); 
for i = 1:N
    if v_i(i) < 0.8                       
       Z_i(i) = 1;  
    else             
       Z_i(i) = 0;
    end
end

X_i = a0 + a1*Z_i + a2*W_i + u_i; 
Y = b0 + b1*X_i + b2*W_i + e_i;

% Primera etapa.

Z = [ones(N,1), Z_i];          % Omitimos W_i.
zeta_gorro = inv(Z'*Z)*Z'*X_i;  

% Segunda etapa.

X_hat = Z*zeta_gorro;

X = [ones(N,1), X_hat];         % Omitimos W_i.
beta_gorro = inv(X'*X)*X'*Y;

% Errores estándares - Segunda Etapa.

Y = b0 + b1*X_i + e_i;
X = [ones(N,1), X_hat];         
beta_gorro_se = inv(X'*X)*X'*Y;
e_gorro_se = Y - X*beta_gorro_se;

K = length(beta_gorro_se);

s = (e_gorro_se'*e_gorro_se)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    % Error estándar homocedástico.
%se_r = sqrt(diag(inv(X'X)(X'(diag(e_gorro_se.^2))(N/(N-K))*X)*inv(X'*X))); % Error estándar robusto.

% Errores estándar - Primera Etapa.

u_gorro = X_i - Z*zeta_gorro;

K = length(zeta_gorro);

s1 = (u_gorro'*u_gorro)/(N-K);    % Varianza estimada del error.
se1 = sqrt(s1*diag(inv(Z'*Z)));    % Error estándar homocedástico.
%se1_r = sqrt(diag(inv(Z'Z)(Z'diag(u_gorro.^2)(N/(N-K))*Z)*inv(Z'*Z))); % Error estándar robusto.

% Test F (Utiliza errores estándares robustos). 

R = [0 1]';
    c = 0;
    q = 1;
    var_beta = se1(2)^2;
    ftest = (R' * zeta_gorro - c)^2 / var_beta;
    p_value1 = 1 - fcdf(ftest, q, N - K);

% Resultados.

Primera_Etapa = table(zeta_gorro, se1, ...
    'VariableNames', {'Coeficiente', 'SE'});

Segunda_Etapa = table(beta_gorro, se, ...
    'VariableNames', {'Coeficiente', 'SE'});

disp(['Para un alpha_1 = ', num2str(a1), ', la primera etapa tiene coeficientes y errores estándares:']);
disp(Primera_Etapa);

disp(['donde el Test F = ', num2str(ftest)]);

disp(['y su segunda etapa será:']);
disp(Segunda_Etapa);

end
%% 4
% Parametros dados en el problema
beta0 = 1;
beta1 = 2;
beta2 = 5;
alpha0 = -4;
alpha2 = 3;
N = 1000;
num_simulations = 1000;
a1_values = [0.1, 0.5, 1, 5, 10];

% Preparar matrices para almacenar resultados
beta1_hat_MC2E = zeros(num_simulations, length(a1_values));
bias_beta1_MC2E = zeros(length(a1_values), 1);

for j = 1:length(a1_values)
    alpha1 = a1_values(j);
    for sim = 1:num_simulations
        % Generación de datos
        W = normrnd(2, 1, N, 1);
        ei = normrnd(0, 1, N, 1);
        ui = normrnd(0, 1, N, 1);
        vi = rand(N, 1);
        Z = double(vi < 0.8);
        
        Xi = alpha0 + alpha1 * Z + alpha2 * W + ui;
        Y = beta0 + beta1 * Xi + beta2 * W + ei;
        
        % Primera etapa
        Z_ext = [ones(N, 1), Z];
        zeta_gorro = (Z_ext' * Z_ext) \ (Z_ext' * Xi);
        X_hat = Z_ext * zeta_gorro;
        
        % Segunda etapa
        X_ext = [ones(N, 1), X_hat];
        b_iv = (X_ext' * X_ext) \ (X_ext' * Y);
        
        % Guardar estimaciones de beta1
        beta1_hat_MC2E(sim, j) = b_iv(2);
    end
    
    % Calcular el sesgo de beta1_hat
    bias_beta1_MC2E(j) = mean(beta1_hat_MC2E(:, j)) - beta1;
    
end



% Graficar la distribución asintótica de beta1_hat para cada alpha1 en un solo gráfico
figure;
hold on;
colors = lines(length(a1_values));
legendEntries = cell(length(a1_values), 1);  % Para almacenar las entradas de la leyenda
for j = 1:length(a1_values)
    h = histogram(beta1_hat_MC2E(:, j), 'Normalization', 'pdf', 'FaceColor', colors(j, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'BinWidth', 0.05);
    % Crear una entrada de leyenda con un rectángulo del color correspondiente
    legendEntries{j} = plot(nan, nan, 's', 'MarkerSize', 10, 'MarkerFaceColor', colors(j, :), 'MarkerEdgeColor', colors(j, :), 'DisplayName', ['\alpha_1 = ', num2str(a1_values(j))]);
end
hold off;
title('Distribución Asintótica de \beta_{1}^{MC2E} para Diferentes Valores de \alpha_{1}');
xlabel('\beta_{1}^{MC2E}');
ylabel('Densidad');
xlim([0 4]); % Ajustar límites de los ejes
legend([legendEntries{:}]);  % Mostrar la leyenda con las entradas personalizadas

% Graficar la distribución asintótica de beta1_hat para cada alpha1 en un solo gráfico


% Reporte del sesgo promedio para cada alpha1
disp('Sesgo promedio de beta1_hat_MC2E para diferentes valores de alpha1:');
for j = 1:length(a1_values)
    disp(['alpha1 = ', num2str(a1_values(j)), ', Bias(beta1_hat_MC2E) = ', num2str(bias_beta1_MC2E(j))]);
end


