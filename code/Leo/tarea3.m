% Tarea 3 - Econometría I
% Nicolás Bastías, Ignacio Leal, Sebastián Mejías y Leonardo Siles

%% 1. Pregunta 1. Probabilidad lineal, logit y probit

% 1.0 Preliminares
clear;
dataset = readtable("data/cps09mar.xlsx");
addpath("code/econometrics/");
y = [dataset.union];
X = [dataset.age dataset.education dataset.hisp dataset.race];
% Muestra de hombres
y = y(dataset.female == 0);
X = X(dataset.female == 0, :);
k = size(X, 2);
% Dummies para cada categoria de raza
% Se hace una revision de la estadística descriptiva para determinar el
% numero de dummies a crear para la variable categórica
tabulate(X(:, 4))
% Donde se tiene que las categorias 1, 2, 3 y 4 agrupan a más del 97.5% de
% toda la muestra. Por tanto, se crean dummies para cada una de estas
% categorías. La dummy para la categoría 'otros' no se crea pues esta sera
% la categoria base de la estimacion
for j = 1:4
    X(:, k+j) = (X(:, k) == j);
end
% Matriz de regresores con intercepto y sin la variable original de raza
X = [X, ones(size(X,1),1)];
X(:, 4) = [];
% Computo de valores de n y k (con las dummies incluidas)
n = size(X,1);
k = size(X,2);
% 1.1 Modelo de probabilidad lineal
% Recordar incrementar el tamaño máximo de arrays en las preferencias de 
% MATLAB
[b_lineal, ~, robust_lineal] = mco(y, X);
serob_lineal = sqrt(diag(robust_lineal));
% 1.2 Probit con el supuesto de especificación correcta
[b_probit, v0_probit, robust_probit, ll_probit] = probit(y, X);
se0_probit = sqrt(diag(v0_probit./n));
% 1.3 Probit con errores estándar robustos
serob_probit = sqrt(diag(robust_probit./n));
% 1.4 Probit - efectos marginales en el promedio
x_means = mean(X);
margins_probit = b_probit.*normpdf(x_means*b_probit);
% Ajuste por variables dicotomicas y categoricas
for j = 3:7
    margins_probit(j) = normcdf(x_means(1:j-1)*b_probit(1:j-1) + ...
        b_probit(j) + x_means(j+1:8)*b_probit(j+1:8)) - normcdf( ...
        x_means(1:j-1)*b_probit(1:j-1) + x_means(j+1:8)*b_probit(j+1:8));
end 
% 1.5 Logit con el supuesto de especificación correcta
[b_logit, v0_logit, robust_logit, ll_logit] = logit(y, X);
se0_logit = sqrt(diag(v0_logit./n));
% 1.6 Logit con errores estándar robustos
serob_logit = sqrt(diag(robust_logit./n));
% 1.7 Logit - efectos marginales en el promedio
pd = makedist('Logistic');
margins_logit = b_logit.*pdf(pd, x_means*b_logit);
% Ajuste por variables dicotomicas y categoricas
for j = 3:7
    margins_logit(j) = cdf(pd, x_means(1:j-1)*b_logit(1:j-1) + ...
        b_logit(j) + x_means(j+1:8)*b_logit(j+1:8)) - cdf(pd, ...
        x_means(1:j-1)*b_logit(1:j-1) + x_means(j+1:8)*b_logit(j+1:8));
end 
% 1.8 Efectos marginales de la educación
edumargins_probit = zeros(size(unique(dataset.education),1), 1);
edumargins_logit = zeros(size(unique(dataset.education),1), 1);
% Note que son efectos marginales evaluados en el promedio, siguiendo el
% enunciado para el ejercicio anterior
for i = unique(dataset.education)
    edumargins_probit = b_probit(2)*normpdf(x_means(1)*b_probit(1) + ...
        i*b_probit(2) +  x_means(3:end)*b_probit(3:end));
     edumargins_logit = b_logit(2)*pdf(pd, x_means(1)*b_logit(1) + ...
        i*b_logit(2) +  x_means(3:end)*b_logit(3:end));
end
% Grafico de efectos marginales del MPL, modelo probit y logit
educ = unique(dataset.education);
edu_mco = ones(size(educ,1),1)*b_lineal(2);
figure;
hold on;
plot(educ, edumargins_probit, '-o', 'LineWidth', 1.5);
plot(educ, edumargins_logit, '-s', 'LineWidth', 1.5);
plot(educ, edu_mco, '-^', 'LineWidth', 1.5);
grid on;
grid minor;
xlabel('Nivel de educación');
ylabel('Efectos marginales');
legend('Probit', 'Logit', 'MPL', 'Location', 'Best');
