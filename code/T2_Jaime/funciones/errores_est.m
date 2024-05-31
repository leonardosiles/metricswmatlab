function [ee,ee_robust,ee_cluster] = errores_est(Y,X,n_groups)

% Estimador de MCO
bhat = (X'*X)^(-1)*X'*Y;

% Residuos de MCO
ehat = Y - X*bhat; 

% Número de regresores y observaciones
N = size(X,1);
K = size(X,2);

% a) Errores estándar homocedasticos
s2 = 1/(N-K)*(ehat'*ehat);
VarB = s2*(X'*X)^(-1);
ee = diag(VarB).^(1/2);

% b) Errores estándar robustos HC1 (corregidos por sesgo de los residuos)
D = diag(ehat.^2);
VarB_robust = N/(N-K)*(X'*X)^(-1)*X'*D*X*(X'*X)^(-1);
ee_robust = diag(VarB_robust).^(1/2);

% c) Errores estándar Agrupados
Omega = zeros(K);
G = n_groups;
size_group = N/G;
for i=1:G
    ind = [i*size_group-(size_group-1) i*size_group];
    X_g = X(ind(1):ind(2),:);
    ehat_g = ehat(ind(1):ind(2));
    aux = X_g'*(ehat_g*ehat_g')*X_g;
    Omega = Omega + aux;
end
VarB_cluster = (G/(G-1))*(X'*X)^(-1)*Omega*(X'*X)^(-1);
ee_cluster = diag(VarB_cluster).^(1/2);

end
