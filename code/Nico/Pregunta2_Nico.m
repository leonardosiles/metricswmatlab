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

X_i = a0 + a1*Z_i + a2*W_i + u_i;    % Primera etapa.
Y = b0 + b1*X_i + a2*W_i + e_i;    % Segunda etapa.

x_i = a0 + a1*Z_i + u_i;    % Primera etapa omitiendo W_i.

X = [ones(N,1), x_i];

beta_gorro = inv(X'*X)*X'*Y;

u_gorro = Y - X*beta_gorro;

K = length(beta_gorro);

% Errores estándares.

s = (u_gorro'*u_gorro)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    % Error estándar homocedástico.
se_r = sqrt(diag(inv(X'*X)*(X'*((diag(u_gorro.^2))/(N-K))*X)*inv(X'*X))); % Error estándar robusto.

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
    Y = b0 + b1*X_i + a2*W_i + e_i;    % Segunda etapa.

    x_i = a0 + a1*Z_i + u_i;    % Primera etapa omitiendo W_i.

    X = [ones(N,1), x_i];
    
    m = randi(N,N,1);      % Vector de m pares aleatorios intependientes de tamaño n. 
    bm(:,i) = (X(m,:)'*X(m,:))\(X(m,:)'*Y(m)) ;
end
bm = sort(bm,2);

% Calcular la densidad estimada
hold on;
[f,xi] = ksdensity(bm(2, :));
plot(xi,f,'LineWidth',2);
hold off;

% Con 1000 muestras aleatorias, el estimador es bien comportado. 

% Sesgo.

sesgo = mean(bm(2,:))- b1;

disp(['Para un alpha_1 = ', num2str(a1), ', el sesgo de nuestro estimador de MCO será:']);
    disp(sesgo);
end

xlabel('Valor');
ylabel('Densidad');
title('Curva de densidad estimada $\hat{\beta}_{1}$', 'Interpreter','latex');
legend('$\alpha_1 = 0.1$', '$\alpha_1 = 0.5$','$\alpha_1 = 1$', ...
    '$\alpha_1 = 5$', '$\alpha_1 = 10$','Interpreter','latex');

filename = ['densidad_alpha', '.png'];
    saveas(gcf, filename);
% El sesgo se reduce con un alpha_1 más grande. 
% Según la el manual de ivregress en stata: Both theoretical and Monte Carlo exercises indicate 
% that the LIML estimator may yield less bias and confidence intervals with better coverage rates 
% than the 2SLS estimator. See Poi (2006) and Stock, Wright, and Yogo (2002) (and the papers cited 
% therein) for Monte Carlo evidence

%% Pregunta 3: Estimación por MC2E.

for a1 = [0.1, 0.5, 1, 5, 10]

% Primera etapa.

X_i = a0 + a1*Z_i + a2*W_i + u_i;    

Z = [ones(N,1), Z_i, W_i];        
zeta_gorro = inv(Z'*Z)*Z'*X_i;

% Test F. Utiliza errores estándares robustos. 

u_gorro = Y - X*beta_gorro;

K = length(zeta_gorro);

se1_r = sqrt(diag(inv(Z'*Z)*(Z'*diag((u_gorro.^2)/(N-K))*Z)*inv(Z'*Z))); % Error estándar robusto (primera etapa)

R = [0 1 0]';
c = [0 0 0]';                      % Valor de la hipotesis a testear -> H0 : a1 = 0.
q = 1;

ftest= (R'*zeta_gorro - c)' * inv(se1_r(2,1)*R'*inv(Z'*Z)*R) * (R'*zeta_gorro - c);  % Teest F.
p_value1 = 2 * (1 - fcdf(abs(ftest), q, N -K)); 

% Segunda etapa.

Y = b0 + b1*X_i + b2*W_i + e_i;    

X = [ones(N,1), X_i, W_i];

beta_gorro = inv(X'*X)*X'*Y;

u_gorro = Y - X*beta_gorro;

K = length(beta_gorro);

% Errores estándares.

s = (u_gorro'*u_gorro)/(N-K);    % Varianza estimada del error.
se = sqrt(s*diag(inv(X'*X)));    % Error estándar homocedástico.
se_r = sqrt(diag(inv(X'*X)*(X'*((diag(u_gorro.^2))/(N-K))*X)*inv(X'*X))); % Error estándar robusto.

Primera_Etapa = table(zeta_gorro, se1_r, ...
    'VariableNames', {'Coeficiente', 'SE Robusto'});

Pregunta3 = table(beta_gorro, se, se_r, ...
    'VariableNames', {'Coeficiente', 'SE', 'SE Robusto'});


disp(['Para un alpha_1 = ', num2str(a1), ', la primera etapa tiene coeficientes y errores estándares:']);
disp(Primera_Etapa);

disp(['donde el Test F = ']);
disp(ftest );

disp(['y su segunda etapa será:']);
disp(Pregunta3);

end

%% Pregunta 4: Distribución asintótica de MC2E.


for a1 = [0.1, 0.5, 1, 5, 10]

a1 = 0.1;

X_i = a0 + a1*Z_i + a2*W_i + u_i;   % Primera etapa.
Y = b0 + b1*X_i + b2*W_i + e_i;    % Segunda etapa.

X = [ones(N,1), X_i, W_i];

B = 1000;         % Número de simulaciones. 
bm = NaN(3,B);    % Vector de coeficientes simulados por Montecarlo.

for i = 1:B
    m = randi(N,N,1);      % Vector de m pares aleatorios intependientes de tamaño n. 
    bm(:,i) = (X(m,:)'*X(m,:))\(X(m,:)'*Y(m)) ;
end
bm = sort(bm,2);

mean0 = mean(bm(1,:))
mean1 = mean(bm(2,:))
mean2 = mean(bm(3,:))

% Calcular la densidad estimada
[f,xi] = ksdensity(bm(2, :));

% Graficar la curva de densidad estimada
plot(xi,f,'LineWidth',2);
xlabel('Valor');
ylabel('Densidad');
title('Curva de densidad estimada $\hat{\beta}_{1}$', 'Interpreter','latex');

filename = ['densidad_alpha_', num2str(a1), '.png'];
    saveas(gcf, filename);

% Con 1000 muestras aleatorias, el estimador es bien comportado. 

% Sesgo.

sesgo = mean(bm(2,:))- b1;

disp(['Para un alpha_1 = ', num2str(a1), ', el sesgo de nuestro estimador de MCO será:']);
    disp(sesgo);
end



