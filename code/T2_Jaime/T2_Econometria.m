%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Tarea Computacional 2                          %
%                              Econometría I                              %
%                      Profesora: Valentina Paredes                       %
%             Ayudantes: Hriday Karnani & María Jesús Negrete             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estudiantes:
% Fernanda Anguita, Jaime Aránguiz, Nicolás Espinoza y Leonardo Siles

%% 0. Importación de base
clear;clc;
cd("C:\Users\jaime\Desktop\FEN\SEMESTRE 9 - OTOÑO 2024 - ME1\Econometría I\Tareas\Tarea2");
addpath("datos");
addpath("funciones");

% importamos la base
txt = importdata("WAGEPAN.txt");
data = txt.data;
variables = txt.colheaders;
n = size(data,1);

% Guardamos los regresores
regresores = {'educ','black','hisp','exper','expersq','married','union'};
X = zeros(n,numel(regresores));
for i=1:numel(regresores)
    index = find(strcmp(regresores(i),variables));
    X(:,i) = data(:,index);
end

% Guardamos la variable dependiente
Y = data(:,strcmp('lwage',variables));

% Guardamos el identificador y los años
ids = data(:,1); N = n/8; T = 8; % tendremos 8 veces la información de n=545 hombres
repeats = diff(ids) == 0;
id = ids;
id(repeats) = [];

% Y los años
year = data(:,2);

%% 1. Pooled MCO
% Estimamos el estimador de pooled MCO
bhat_pool = (X'*X)^(-1)*X'*Y;
display(bhat_pool')

% Guardamos los residuos y la cantidad de regresores
ehat_pool = Y - X*bhat_pool; 
K = size(X,2);

% Calculamos cada error estándar con la función que creamos previamente
[ee_pool,ee_r_pool,ee_cl_pool] = errores_est(Y,X,N);

%% 2. Estimación Between
% Calculamos los promedios a través del tiempo a nivel individual
Y_mean = zeros(N,1);
X_mean = zeros(N,K);
for i=1:N
    Y_mean(i) = mean(Y((T*i-(T-1)):(T*i)));
    X_mean(i,:) = mean(X((T*i-(T-1)):(T*i),:),1);
end

% Estimamos el estimador between
bhat_bet = (X_mean'*X_mean)^(-1)*X_mean'*Y_mean;
display(bhat_bet')

% Guardamos los residuos
ehat_bet = Y_mean - X_mean*bhat_bet; 

% Calculamos cada error estándar con la función que creamos previamente
[ee_bet,ee_r_bet,ee_cl_bet] = errores_est(Y_mean,X_mean,N);

%% 4. Estimación de efectos fijos (within)
% Realizaremos primero la estimación mediante efectos fijos para luego
% poder construir el estimador de efectos aleatorios

% Eliminamos los regresores que no varían en el tiempo
X_fe = X(:,4:end);

% Hacemos la transformación within
Y_tilde = zeros(n,1);
X_tilde = zeros(n,size(X_fe,2));
aux_meanY = repelem(Y_mean,T);
aux_meanX = repelem(X_mean(:,4:end),T,1);
for i=1:n
    Y_tilde(i) = Y(i) - aux_meanY(i);
    X_tilde(i,:) = X_fe(i,:) - aux_meanX(i,:);
end

bhat_fe = (X_tilde'*X_tilde)^(-1)*X_tilde'*Y_tilde;
display(bhat_fe')

% Forma alternativa
%{
% Generamos las dummies de efecto fijo individual
dummies = zeros(n,N);
for i=1:length(id)
    aux = id(i);
    dummies(ids == aux,i) = 1;
end

% Estimamos el estimador de efectos fijos mediante regresión particionada
% Creamos la matriz de aniquilación
M = eye(size(dummies,1)) - dummies*(dummies'*dummies)^(-1)*dummies';

% Eliminamos los regresores que no varían en el tiempo
X_fe = X(:,4:end);

% Estimamos el estimador de efectos fijos
bhat_fe = (X_fe'*M*X_fe)^(-1)*X_fe'*M*Y;
display(bhat_fe)
%}

% Guardamos los residuos
ehat_fe = Y_tilde - X_tilde*bhat_fe; 

% Calculamos cada error estándar con la función que creamos previamente
[ee_fe,ee_r_fe,ee_cl_fe] = errores_est(Y_tilde,X_tilde,N);

%% 3. Estimación de efectos aleatorios
% Ahora podemos estimar mediante efectos aleatorios
% Primero calculamos los parámetros sigma_u y sigma_u
% Sigma_e
suma = 0;
for i=1:n
    suma = suma + ehat_fe(i)^2;
end
sigma_e = 1/(n-N-K)*suma;

% Sigma_u
suma = 0;
for i=1:N
    suma = suma + ehat_bet(i)^2;
end
sigma_u = 1/(N-K)*suma - 1/T *sigma_e;

% Calculamos el rho para poder estimar el Mínimos Cuadrados Generalizados
% Factibles (Feasible GLS)
rho = sigma_e^(1/2)/sqrt(sigma_e + T * sigma_u);

% Estimamos mediante Mínimos Cuadrados Generalizados Factibles
% Para esto, hacemos la transformación de la sección 17.15 de Hansen (2022)
Y_tilde = zeros(n,1);
X_tilde = zeros(n,size(X,2));
aux_meanY = repelem(Y_mean,T);
aux_meanX = repelem(X_mean,T,1);
for i=1:n
    Y_tilde(i) = Y(i) - (1-rho)*aux_meanY(i);
    X_tilde(i,:) = X(i,:) - (1-rho)*aux_meanX(i,:);
end

% Estimamos el estimador de efectos aleatorios
bhat_re = (X_tilde'*X_tilde)^(-1)*X_tilde'*Y_tilde;
display(bhat_re')

% Guardamos los residuos
ehat_re = Y_tilde - X_tilde*bhat_re; 

% Calculamos cada error estándar con la función que creamos previamente
[ee_re,ee_r_re,ee_cl_re] = errores_est(Y_tilde,X_tilde,N);

%% 5. Intervalos de confianza
% Extraemos el error estándar y el estimador para el premio de pertenecer a
% un sindicato

%%%%%%%%%%%%%%%%%%%%%%%%%%   TEORÍA ASINTÓTICA   %%%%%%%%%%%%%%%%%%%%%%%%%%
% Guardamos el estimador y error estándar de beta7
se = ee_cl_pool(end);
bhat = bhat_pool(end);

% Calculamos el intervalo de cofianza al 95% asumiendo que el estimador se
% distribye normal
ci_asintotica = [(bhat-1.96*se) (bhat+1.96*se)];

%%%%%%%%%%%%%%%%%%%%%%%%   BLOCK BOOTSTRAP   %%%%%%%%%%%%%%%%%%%%%%%%%
rng(7);                     % Fijamos una seed
n_iter = 1000;              % Fijamos un número de replicaciones
b_block = zeros(n_iter,1);   % Inicializamos la matriz donde vamos a guardar los estimadores de cada replicación

% Efectuamos bootstrap
tic
for r=1:n_iter
    % Generamos las observaciones de bootstap muestreando individuos
    index = randsample(id,N,true);
    
    Y_w = zeros(n,1); 
    X_w = zeros(n,K);
    ehat_w = zeros(n,1);
    for i=1:N
        Y_w((T*i-(T-1)):(T*i)) = Y(ids == index(i));
        X_w((T*i-(T-1)):(T*i),:) = X(ids == index(i),:);
        ehat_w((T*i-(T-1)):(T*i)) = ehat_pool(ids == index(i));
    end

    % Calculamos el estimador de bootstrap
    bhat_wild = (X_w'*X_w)^(-1)*X_w'*Y_w;

    % Almacenamos beta7 en la matriz de estimadores
    b_block(r) = bhat_wild(end);
end
toc

% Calculamos el intervalo de confianza mediante el método del percentil
ci_blockboot = [prctile(b_block,2.5) prctile(b_block,97.5)];

%%%%%%%%%%%%%%%%%%%%%%%   WILD CLUSTER BOOTSTRAP   %%%%%%%%%%%%%%%%%%%%%%%%

% Inicializamos la matriz donde vamos a guardar los estimadores de cada replicación
b_cluster = zeros(n_iter,1);

% Efectuamos bootstrap
tic
for r=1:n_iter
    index = randsample(id,N,true);

    Y_w = zeros(n,1); 
    X_w = zeros(n,K);
    ehat_w = zeros(n,1);
    for i=1:N
        Y_w((T*i-(T-1)):(T*i)) = Y(ids == index(i));
        X_w((T*i-(T-1)):(T*i),:) = X(ids == index(i),:);
        ehat_w((T*i-(T-1)):(T*i)) = ehat_pool(ids == index(i));
    end
    
    % Generamos la variable auxiliar xi con la distribución propuesta por
    % Mammen (1993)
    p = binornd(1,(sqrt(5)-1)/(2*sqrt(5)),[N 1]);
    xi = p*(1+sqrt(5))/2 + (1-p)*(1-sqrt(5))/2;
    
    % Ajustamos el residuo por la variable auxiliar
    ehat_star = zeros(n,1);
    for i=1:N
        ehat_star((T*i-(T-1)):(T*i)) = ehat_w((T*i-(T-1)):(T*i)).*xi(i);
    end
    
    % Con el residuo ajustado generamos la observación de bootstrap
    % definitiva para Y
    Y_star = X_w*bhat_pool + ehat_star; 
    
    % Calculamos el estimador de bootstrap
    bhat_wild = (X_w'*X_w)^(-1)*X_w'*Y_star;
    
    % Almacenamos beta7 en la matriz de estimadores
    b_cluster(r) = bhat_wild(end);
end
toc

% Calculamos el intervalo de confianza mediante el método del percentil
ci_clusterboot = [prctile(b_cluster,2.5) prctile(b_cluster,97.5)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RESULTADOS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table(ci_asintotica,ci_blockboot,ci_clusterboot)

%% 6. Incluyendo efectos fijos temporales - Pooled MCO
% Guardamos las dummies de cada año en una nueva matriz "v"
dummies_year = {'d81','d82','d83','d84','d85','d86','d87'};
v = zeros(n,numel(dummies_year));
for i=1:numel(dummies_year)
    index = find(strcmp(dummies_year(i),variables));
    v(:,i) = data(:,index);
end

% Definimos una nueva matriz de regresores con los efectos fijos anuales
Xv = [X v];

% Estimador between
bhat2_pool = (Xv'*Xv)^(-1)*Xv'*Y;
display(bhat2_pool(1:7)');

% Calculamos cada error estándar con la función que creamos previamente
[ee2_pool,ee2_r_pool,ee2_cl_pool] = errores_est(Y,Xv,N);


%% 7. Incluyendo efectos fijos temporales - Fixed effects
% Calculamos promedios individuales para cada regresor
Xv_mean = zeros(N,size(Xv,2));
for i=1:N
    Xv_mean(i,:) = mean(Xv((T*i-(T-1)):(T*i),:),1);
end

% Eliminamos los regresores que no varían en el tiempo
Xv_fe = Xv(:,4:end);

% Hacemos la transformación within
Y_tilde = zeros(n,1);
Xv_tilde = zeros(n,size(Xv_fe,2));
aux_meanY = repelem(Y_mean,T);
aux_meanX = repelem(Xv_mean(:,4:end),T,1);
for i=1:n
    Y_tilde(i) = Y(i) - aux_meanY(i);
    Xv_tilde(i,:) = Xv_fe(i,:) - aux_meanX(i,:);
end

%bhat2_fe = (Xv_tilde'*Xv_tilde)^(-1)*Xv_tilde'*Y_tilde;
[bhat2_fe, ~, y_ddot, x_ddot] = twfe(Y, X, ids, year);
display(bhat2_fe')


% Calculamos cada error estándar con la función que creamos previamente
[ee2_fe,ee2_r_fe,ee2_cl_fe] = errores_est(Y_tilde,Xv_tilde,N);

%% 8. Testeando cambios en el retorno de la educación
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ESTIMACIÓN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primero creamos las interacciones entre las dummies anuales y la variable
% de educación
interacciones = zeros(n,size(v,2));
for k=1:size(v,2)
    interacciones(:,k) = X(:,1).*v(:,k);
end

% Definimos la nueva matriz de regresores, incluyendo interacciones
Xv_int = [Xv interacciones];

% Calculamos promedios individuales para cada regresor
Xv_int_mean = zeros(N,size(Xv_int,2));
for i=1:N
    Xv_int_mean(i,:) = mean(Xv_int((T*i-(T-1)):(T*i),:),1);
end

% Eliminamos los regresores que no varían en el tiempo
Xv_int_fe = Xv_int(:,4:end);

% Hacemos la transformación within para el regresor (para Y ya la tenemos)
Xv_int_tilde = zeros(n,size(Xv_int_fe,2));
aux_meanX = repelem(Xv_int_mean(:,4:end),T,1);
for i=1:n
    Xv_int_tilde(i,:) = Xv_int_fe(i,:) - aux_meanX(i,:);
end

% Calculamos el estimador within
bhat_interac = (Xv_int_tilde'*Xv_int_tilde)^(-1)*Xv_int_tilde'*Y_tilde;

% Guardamos los residuos
ehat_interac = Y_tilde - Xv_int_tilde*bhat_interac; 

% Calculamos la matriz de varianza y covarianzas de datos agrupados a mano,
% para obtener la matriz de varianza y covarianzas del estimador
K = size(Xv_int_tilde,2);
Omega = zeros(K);
G = N;
size_group = n/G;
for i=1:G
    ind = [i*size_group-(size_group-1) i*size_group];
    X_g = Xv_int_tilde(ind(1):ind(2),:);
    ehat_g = ehat_interac(ind(1):ind(2));
    aux = X_g'*(ehat_g*ehat_g')*X_g;
    Omega = Omega + aux;
end
Var_b = (G/(G-1))*(n-1)/(n-K)*(Xv_int_tilde'*Xv_int_tilde)^(-1)*Omega*(Xv_int_tilde'*Xv_int_tilde)^(-1);
%ee_interac_cl = diag(VarB_cluster).^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%   TEST DE HIPÓTESIS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definimos el q del test conjunto
q = 7;

% Repetimos Definimos el vector c
c = repelem(0,q,1);

% Calculamos R'b
R = [zeros(11,7); eye(7)];
Rb = R'*bhat_interac;

% Calculamos el estadístico Generalizado de Wald del test conjunto
%F_statistic = (Rb-c)'*(R'*Var_b*R)^(-1)*(Rb-c) /q;
GW_statistic = (Rb-c)'*(R'*Var_b*R)^(-1)*(Rb-c);

% Calculamos el p-value correspondiente
%p_value = 1 - fcdf(F_statistic,q,n-K);
p_value = 1 - chi2cdf(GW_statistic,q);

% Mostramos el resultado
table(GW_statistic, p_value)

%% 9. Testeando cambios en el premio por sindicato según raza
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ESTIMACIÓN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primero creamos las interacciones entre las variables de raza (black e
% hisp) y la variable de sindicato (union)
interacciones = zeros(n,2);             % inicializamos la matriz de interacciones
interacciones(:,1) = X(:,7).*X(:,2);    % Creamos union * black
interacciones(:,2) = X(:,7).*X(:,3);    % Creamos union * hisp

% Definimos la nueva matriz de regresores, incluyendo interacciones
Xv_int = [Xv interacciones];

% Calculamos promedios individuales para cada regresor
Xv_int_mean = zeros(N,size(Xv_int,2));
for i=1:N
    Xv_int_mean(i,:) = mean(Xv_int((T*i-(T-1)):(T*i),:),1);
end

% Eliminamos los regresores que no varían en el tiempo
Xv_int_fe = Xv_int(:,4:end);

% Hacemos la transformación within para el regresor (para Y ya la tenemos)
Xv_int_tilde = zeros(n,size(Xv_int_fe,2));
aux_meanX = repelem(Xv_int_mean(:,4:end),T,1);
for i=1:n
    Xv_int_tilde(i,:) = Xv_int_fe(i,:) - aux_meanX(i,:);
end

% Calculamos el estimador within
bhat_interac2 = (Xv_int_tilde'*Xv_int_tilde)^(-1)*Xv_int_tilde'*Y_tilde;

% Guardamos los residuos
ehat_interac2 = Y_tilde - Xv_int_tilde*bhat_interac2; 

% Calculamos la matriz de varianza y covarianzas de datos agrupados a mano,
% para obtener la matriz de varianza y covarianzas del estimador
K = size(Xv_int_tilde,2);
Omega = zeros(K);
G = N;
size_group = n/G;
for i=1:G
    ind = [i*size_group-(size_group-1) i*size_group];
    X_g = Xv_int_tilde(ind(1):ind(2),:);
    ehat_g = ehat_interac2(ind(1):ind(2));
    aux = X_g'*(ehat_g*ehat_g')*X_g;
    Omega = Omega + aux;
end
Var_b = (G/(G-1))*(N-1)/(N-K)*(Xv_int_tilde'*Xv_int_tilde)^(-1)*Omega*(Xv_int_tilde'*Xv_int_tilde)^(-1);
%ee_interac_cl = diag(VarB_cluster).^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%   TEST DE HIPÓTESIS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testearemos la hipótesis nula de que los b de las interacciones son = 0
% Definimos el q del test conjunto
q = 2;

% Creamos el vector c = (0 0)'
c = repelem(0,q,1);

% Calculamos R'b
R = [zeros(11,2); eye(2)];
Rb = R'*bhat_interac2;

% Calculamos el estadístico Generalizado de Wald del test conjunto
GW_statistic = (Rb-c)'*(R'*Var_b*R)^(-1)*(Rb-c);

% Calculamos el p-value correspondiente
p_value = 1 - chi2cdf(GW_statistic,q);

% Mostramos el resultado
table(GW_statistic, p_value)

%% 10. Agregando un adelanto a la ecuación 2
% Primero debemos la nueva matriz de regresores que incluya el vector con 
% el adelanto de la variable union para cada individuo (excluyendo el 
% último periodo)
K = size(Xv,2);
T2 = T-1;   % Ahora el panel pasará a tener 7 periodos por individuo
n2 = n-N;   % Y n-N = 3815 observaciones totales
Xv_adel = zeros(n2,K);    % Inicializamos la nueva matriz de regresores que K variables (+1-1)
Y_adel = zeros(n2,1);       % Inicializamos el nuevo vector de Y
for i=1:N
    % Extraemos la información de la persona i
    X_i = Xv((T*i-(T-1)):(T*i),:);
    Y_i = Y((T*i-(T-1)):(T*i),1);
    
    % Guardamos sus primeras 7 observaciones de cada periodo para X e Y
    Y_adel((T2*i-(T2-1)):(T2*i),1) = Y_i(1:7,1);
    Xv_adel((T2*i-(T2-1)):(T2*i),1:(K-1)) = X_i(1:7,1:(K-1));   % EXCLUIMOS LA ÚLTIMA DUMMY DE AÑO

    % Guardamos en la última columna el adelanto de union (variable 7)
    Xv_adel((T2*i-(T2-1)):(T2*i),end) = X_i(2:8,7);
end

% Calculamos los promedios de las matrices X e Y
Xv_adel_mean = zeros(N,size(Xv_adel,2));
Y_adel_mean = zeros(N,1);
for i=1:N
    Xv_adel_mean(i,:) = mean(Xv_adel((T2*i-(T2-1)):(T2*i),:),1);
    Y_adel_mean(i,1) = mean(Y_adel((T2*i-(T2-1)):(T2*i),1),1);
end

% Eliminamos los regresores que no varían en el tiempo
Xv_adel_fe = Xv_adel(:,4:end);

% Hacemos la transformación within
Y_tilde = zeros(n2,1);
Xv_tilde = zeros(n2,size(Xv_adel_fe,2));
aux_meanY = repelem(Y_adel_mean,T2);
aux_meanX = repelem(Xv_adel_mean(:,4:end),T2,1);
for i=1:n2
    Y_tilde(i) = Y_adel(i) - aux_meanY(i);
    Xv_tilde(i,:) = Xv_adel_fe(i,:) - aux_meanX(i,:);
end

% Calculamos el estimador within
bhat2_adel = (Xv_tilde'*Xv_tilde)^(-1)*Xv_tilde'*Y_tilde;
display([bhat2_adel(1:4)' bhat2_adel(end)])

% Calculamos cada error estándar con la función que creamos previamente
[ee2_adel,ee2_adel_r,ee2_adel_cl] = errores_est(Y_tilde,Xv_tilde,N);


