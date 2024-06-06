%% Pregunta 2: Tarea 3 (Econometría I)

clear all
clc

rng(14)               % Semilla. 

b0 = 1;
b1 = 2;
b2 = 5;

a0 = -4;
a1 = 0.1;                   % a1 = {0.1, 0.5, 1, 5, 10}
a2 = 3; 

N = 1000;

e_i = 1*randn(N,1);
u_i = 1*randn(N,1);

W_i = 1*randn(N,1) + 2;

v_i = unifrnd(0,1,[N,1]); 

Z_i = zeros(N,1); 

for i = 1:N
    if v_i(i) < 0.8                       
       Z_i(i) = 1;  
    else             
       Z_i(i) = 0;
    end
end

%% Pregunta 1: W_i no observable.

% Si W_i no es observable, tampoco entrará dentro de nuestro instrumento.
% De esta forma, la primera y segunda etapa serán respectivamente:

X_i = a0 + a1*Z_i + u_i;    % Primera etapa.
Y_i = b0 + b1*X_i + e_i;    % Segunda etapa.

X = [ones(N,1), X_i];

beta_gorro = inv(X'*X)*X'*Y_i

e_gorro = Y_i - X*beta_gorro;

    K = length(beta_gorro)

% Errores estándares.

s = (e_gorro'*e_gorro)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    

Pregunta1 = [beta_gorro, se]

% Hacer con un loop para cada alpha. 

