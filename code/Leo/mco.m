% Estimador de MCO
% Incluye tres tipos de errores estandar:
% Homocedasticos, robustos y ajustados por clúster

function [coef, vcov, rvcov, clvcov] = mco(y, X, id)
arguments
    y {mustBeNumeric}
    X {mustBeNumeric}
    id = 0
end
    n = size(X, 1);         % cantidad de observaciones
    k = size(X, 2);         % numero de regresores
    XX = X'*X;              % multiplicacion x transpuesta por x
    %IXX = inv(XX);         % inversa de x transpuesta por x
    coef = XX\(X'*y);       % estimadores de MCO
    e = y - X*coef;         % residuos de MCO

    % Hasta aquí, se calcularon los betas de MCO. Ahora escribimos el
    % código para calcular los tres diferentes tipos de errores estándar
    
    % Errores estándar clásicos 

    ee = e'*e;           % suma de los residuos de MCO al cuadrado, escalar
    sd = ee/(n - k);     % estimador insesgado de la varianza de los errores
    vcov = sd*eye(k)/XX; % matriz de varianzas y covarianzas de coefs MCO


    % Errores estándar robustos a la heterocedasticidad

    e2 = e.^2;              % residuos de MCO al cuadrado, vector n x 1
    %D = diag(e2);           % matriz diagonal con residuos al cuadrado, n x n
    sdw = XX\X'*diag(e2)*X/XX;     % sandwich de White
    rvcov = (n/(n-k))*sdw;  % ajuste por grados de libertad. matriz de varianzas y covarianzas robusta

    % Errores estandar ajustados por cluster
    if id ~= 0
    [~, ~, index] = unique(id);
    g = length(unique(id));
    sumxe = zeros(g,k);     
    for i = 1:k
        sumxe(:,i) = accumarray(index, X(:, i).*e)';
    end
    omega = sumxe'*sumxe;                       % Estimador varianza error
    adj = (g/(g-1))*((n-1)/(n-k));              % Ponderador 
    clvcov = adj*eye(k)*XX\omega/XX;            % Matriz var-cov agrupada 
    end
end

