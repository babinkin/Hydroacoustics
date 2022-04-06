clc
close all
clear 
%==========================================================================
% Настройка шрифтов в зависимости от ОС
if isunix
    fontname = 'Free Helvetian';% for LINUX
    % Для Юникс системвыставляется не тот шрифт в функции msgbox, поэтому
    % либо правим файл toolbox/matlab/uitools/msgbox.m, либо для каждого
    % объекта после создания меняем шрифт
elseif ispc
    fontname = 'Arian Cyr';% for Windows
end
set(0,'DefaultAxesFontName',fontname);
set(0,'DefaultTextFontName',fontname);
set(0,'DefaultUIControlFontname',fontname);
% set(0,'fixedwidthfontname',fontname);
TextSize = 20;
%==========================================================================


logsDir = 'foam_data/logs';
addpath(logsDir)

% -------------------------------------------------------------------------
% Считываем количество итераций на каждом временном шаге
nPimpleIter = dlmread('nPimpleIter'); 
% Количество завершенных шагов по премени
nTimeStep = length(nPimpleIter);

% Считываем текущее время
time_all = dlmread('Time');
time = time_all(1:nTimeStep);

% -------------------------------------------------------------------------
% Считываем данные для давления
p_all = dlmread('p');
p_InitialRes_all = p_all(:,1);
p_FinalRes_all = p_all(:,2);
p_NoIter_all = p_all(:,3);

% Количество итерации давления по неортогональности
% nPressureSolve = nNonOrthogonalCorrectors +1
nPressureSolve = 3;

% Выборка давления
for i=1:nPressureSolve
    p_InitialRes(:,i) = p_InitialRes_all((0 + i):nPressureSolve:(nTimeStep*nPressureSolve));
    p_FinalRes(:,i) = p_FinalRes_all((0 + i):nPressureSolve:(nTimeStep*nPressureSolve));
    p_NoIter(:,i) = p_NoIter_all((0 + i):nPressureSolve:(nTimeStep*nPressureSolve));
end

% -------------------------------------------------------------------------
% Считываем данные для скорости
Ux_all = dlmread('Ux');
Ux_InitialRes = Ux_all(1:nTimeStep,1);
Ux_FinalRes = Ux_all(1:nTimeStep,2);
Ux_NoIter = Ux_all(1:nTimeStep,3);

% --
Uy_all = dlmread('Uy');
Uy_InitialRes = Uy_all(1:nTimeStep,1);
Uy_FinalRes = Uy_all(1:nTimeStep,2);
Uy_NoIter = Uy_all(1:nTimeStep,3);

% --
Uz_all = dlmread('Uz');
Uz_InitialRes = Uz_all(1:nTimeStep,1);
Uz_FinalRes = Uz_all(1:nTimeStep,2);
Uz_NoIter = Uz_all(1:nTimeStep,3);

% -------------------------------------------------------------------------
% Считываем данные для турбулентности
k_all = dlmread('k');
k_InitialRes = k_all(1:nTimeStep,1);
k_FinalRes = k_all(1:nTimeStep,2);
k_NoIter = k_all(1:nTimeStep,3);

% --
omega_all = dlmread('omega');
omega_InitialRes = omega_all(1:nTimeStep,1);
omega_FinalRes = omega_all(1:nTimeStep,2);
omega_NoIter = omega_all(1:nTimeStep,3);




%%

semilogy(1:nTimeStep, p_InitialRes(:,nPressureSolve) ,'-k',1:nTimeStep, Ux_InitialRes ,'-r', 1:nTimeStep, Uy_InitialRes ,'-g', 1:nTimeStep, Uz_InitialRes ,'-b',...
    1:nTimeStep, k_InitialRes ,'-r+',1:nTimeStep, omega_InitialRes ,'-b+','Linewidth',1)
% legend({'p','U_x','U_y','U_z','k','\omega'},'NumColumns',2)
xlim([(nTimeStep - 20) nTimeStep])
% grid on
% xlabel('\alpha, deg','FontSize',20)
set(gca,'FontSize',20)

























