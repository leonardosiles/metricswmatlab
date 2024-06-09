%% Pregunta 2: Tarea 3 (Econometría I)

clear all
clc

rng(14)               % Semilla. 

b0 = 1;
b1 = 2;
b2 = 5;

a0 = -4;
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

for a1 = [0.1, 0.5, 1, 5, 10]

X_i = a0 + a1*Z_i + u_i;    % Primera etapa.
Y = b0 + b1*X_i + e_i;    % Segunda etapa.

X = [ones(N,1), X_i];

beta_gorro = inv(X'*X)*X'*Y;

e_gorro = Y - X*beta_gorro;

K = length(beta_gorro);

% Errores estándares.

s = (e_gorro'*e_gorro)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    

Pregunta1 = [beta_gorro, se];

disp(['Para un \alpha_1 = ', num2str(a1), ', los coeficientes y errores estándares son:']);
    disp(Pregunta1);

end


%% Pregunta 2: Distribución asintótica por Montecarlo.

% Se repite el mismo procedimiento de la pregunta a), pero considerando las
% B = 1000 simulaciones para beta.

B = 1000;         % Número de simulaciones. 
bm = NaN(2,B);    % Vector de coeficientes simulados por Montecarlo.

for i = 1:B
    m = randi(N,N,1);      % Vector de m pares aleatorios intependientes de tamaño n. 
    bm(:,i) = (X(m,:)'*X(m,:))\(X(m,:)'*Y(m)) ;
end
bm = sort(bm,2);

% Calcular la densidad estimada
[f,xi] = ksdensity(bm(2, :));

% Graficar la curva de densidad estimada
plot(xi,f,'LineWidth',2);
xlabel('Valor');
ylabel('Densidad');
title('Curva de densidad estimada $\hat{\beta}_{1}$', 'Interpreter','latex');

% Con 1000 muestras aleatorias, el estimador es bien comportado. 

% Sesgo.

sesgo = mean(bm(2,:))- b1;

disp('El sesgo de nuestro estimador de MCO será:');
    disp(sesgo);

%% Pregunta 3: Estimación por MC2E.





