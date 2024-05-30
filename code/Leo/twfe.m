% Two way fixed effects estimator for balanced panel data
% v. 0.1

function [beta, res, y_ddot, x_ddot] = twfe(Y, X, id, t) 

% Preliminares
N = size(X, 1);
K = size(X, 2);
T = length(unique(t));
n = length(unique(id));
% Primero, asigno indices para los individuos y los años. accumarray
% funciona con esos indices
[~, ~, id] = unique(id);
[~, ~, t] = unique(t);
ym_i = accumarray(id, Y, [], @mean);
ym_t = accumarray(t, Y, [], @mean);
ym_o = mean(Y);
% Inicializo matrices para las medias de la matriz de regresores x
xm_i = zeros(n, K);
xm_t = zeros(T, K); 
xm_o = zeros(1, K);
% Ahora, un for loop para llenar cada columna de las matrices de medias de
% x
for j = 1:K
    xm_i(:, j) = accumarray(id, X(:, j), [], @mean);
    xm_t(:, j) = accumarray(t, X(:, j), [], @mean);
    xm_o(j)    = mean(X(:, j));
end
% Una vez se tienen calculados los tres tipos de medias, se hace la
% transformacion within para TWFE:
y_ddot = Y - ym_i(id) - ym_t(t) + ym_o;
x_ddot = zeros(size(X));
for j = 1:K
    x_ddot(:,j) = X(:,j) - xm_i(id, j) - xm_t(t, j) + xm_o(j); 
end
% Antes de calcular los betas, x_ddot debe ser invertible. El siguiente
% for loop y condicional if verifican esto
indices = zeros(1: K);
for j = 1:K
    temp = sum(x_ddot(:, j));
    if temp == 0
        fprintf("La columna %d de X solo contiene ceros. Por tanto, será eliminada de la matriz\n", j); 
    else 
        indices(j) = j;
    end
end
x_ddot = x_ddot(:, nonzeros(indices)');

% Finalmente, se calculan los coefs y residuos del estimador TWFE
xx = x_ddot'*x_ddot;
beta = xx\x_ddot'*y_ddot;
res = y_ddot - x_ddot*beta;
