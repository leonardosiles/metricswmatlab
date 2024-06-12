% Estimador de modelos probit
% Utiliza el algoritmo de BHHH
% v 0.1

function [b, v0, robust] = probit(y, X, b0, tol)
% Asigno valores por defecto para el beta inicial y el criterio de
% tolerancia
arguments
    y {mustBeNumericOrLogical}
    X {mustBeNumeric}
    b0 = (X'*X)\X'*y
    tol = 1*10^(-9)
end
maxiter = 10^(5);
n = size(X, 1);
k = size(X, 2);
% Ajusto las filas de la matriz X según el valor de la dicotómica y. Guardo
% los resultados de la transformación en una nueva matriz Z
Z = zeros(n, k);
for i = 1:n
    if y(i) == 0
        Z(i, :) = (-1).*X(i, :);
    else
        Z(i, :) = X(i, :);
    end
end
% Inicializo matrices y vectores necesarios para el algoritmo de BHHH
scores = zeros(n, k);
%ng = ones(k, 1);
% Matriz del producto exterior de los scores
B = zeros(k, k);
% Criterio de convergencia. Lo inicializamos en 1
m = 1;

iter = 0;
tic;
while m > tol
    for i = 1:n
        % Scores individuales
        scores(i, :) = (normpdf(Z(i, :)*b0)/normcdf(Z(i, :)*b0)).*Z(i, :);
        % Producto exterior de los scores individuales kxk
        ssi = scores(i,:)'*scores(i,:);
        B = B + ssi;
    end
    % Suma de scores individuales para obtener el gradiente
    g = sum(scores)';
    % Actualización del valor de los betas según el algoritmo BHHH
    b1 = b0 + B\g;
    % Log-verosimilitudes asociadas a b0 y b1
    ll0 = sum(log(normcdf(Z*b0)));
    ll1 = sum(log(normcdf(Z*b1)));
    % Si la log-verosimilitud aumenta con b1, se utiliza el tamaño de salto
    % lambda:
    lambda = 1;
    while ll1 > ll0
         % Duplico el valor de lambda
         lambda = lambda*2;
         % Calcula el nuevo beta asociado al valor duplicado de lambda
         b_lambda1 = b0 + lambda.*(B\g);
         % Calcula las log-verosimilitudes del lambda anterior y actual
         ll0 = ll1;
         ll1 = sum(log(normcdf(Z*b_lambda1)));
    end
    % Finalmente, se calcula b1 con el lambda que hace el mejor salto
    b1 = b0 + lambda.*(B\g);
    % Criterio de convergencia. Estadistico m
    m = g'*(B\g);
    % Actualización de b0 para la siguiente iteración
    b0 = b1;
    %ng = g./b1;
    iter = iter + 1;
    if iter > maxiter
        fprintf('Se ha alcanzado el número máximo de iteraciones (%d). No fue posible lograr la convergencia', maxiter)
        break
    end
    if mod(iter, 5) == 0
        fprintf('El número de iteraciones es %d.\n', iter)
    end
end
toc;
% Beta de probit
b = b1;

%% Queda pendiente el calculo de las matrices de covarianza
v0 = zeros(k, k);
robust = zeros(k, k);
end
