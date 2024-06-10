%% Pregunta 1: Tarea 3 (Econometría I).

datos = xlsread("cps09mar.xlsx");

    age = datos(:,1);
    female = datos(:,2);
    hisp = datos(:,3);
    education = datos(:,4);
    earnings = datos(:,5);
    hours = datos(:,6);
    week = datos(:,7);
    union = datos(:,8);
    uncov = datos(:,9);
    region = datos(:,10);
    race = datos(:,11);
    marital = datos(:,12);
