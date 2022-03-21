clear
close all
clc

% загрузка файла данных

% из Fluent
% data_file = 'point_x_0_1_y_0_025_z_0.dat';
% title_line = 4; % количество строк заголовка

% из FOAM
data_file = 'p';
title_line = 20; % количество строк заголовка

% -------------------------------------------------------------------------
% Количество строк в файле
fileID = fopen(data_file);
Cell_scan = textscan(fileID,'%s','Delimiter','\n');
num_line = length(Cell_scan{1,1});

fclose(fileID);
clear Cell_scan

% -------------------------------------------------------------------------
% считываем данные, без заголовка и последней строки
fileID = fopen(data_file);

% из Fluent
% p_over_t = dlmread(data_file,' ',[title_line 0 (num_line - title_line + 2) 1]);

% из Foam
p_over_t_cell = textscan(fileID,repmat('%f',1,19), 'MultipleDelimsAsOne',true, 'HeaderLines',title_line);
p_over_t(:,1) = p_over_t_cell{1};
p_over_t(:,2) = p_over_t_cell{3}*998.2;

fclose(fileID);

% -------------------------------------------------------------------------
t = p_over_t(:,1) - p_over_t(1,1); % отсчет от 0
dt = t(2) - t(1);


L = length(p_over_t(:,1)); % number sample
N = floor(L/2); % количество частот четное

% -------------------------------------------------------------------------
% Фильтрованая функция

fs = 1/(dt); % data sampled frequency
fc = 300; % cutoff frequency

% Lowpass Butterworth  filter with a cutoff frequency of fc Hz, which, for data sampled at fs Hz
% [b,a] = butter(10,fc/(fs/2),'low');
% [b,a] = cheby2(6,40,fc/(fs/2),'low');

% bandstop filter
% [b,a] = butter(3,[50/(fs/2) 250/(fs/2)],'stop'); 

% highpass Butterworth filter
% [b,a] = butter(5,10/(fs/2),'high');
[b,a] = cheby2(6,60,60/(fs/2),'high');


% Butterworth bandpass filter
% [b,a] = butter(3,[100/(fs/2) 200/(fs/2)],'bandpass'); 
% [b,a] = cheby2(4,40,[100/(fs/2) 200/(fs/2)],'bandpass');

figure(3)
freqz(b,a)
% p_over_t(:,2) = filter(b,a,p_over_t(:,2));



% -------------------------------------------------------------------------
p_mean = mean(p_over_t(:,2));
p_min = min(p_over_t(:,2));
p_max = max(p_over_t(:,2));
p = p_over_t(:,2) - p_mean;
p_rms = rms(p);

f_min_data = 1/(L*dt) %  Минимально возможная частота из сигнала по 2м точнкам

% -------------------------------------------
% график давления от 
figure(1)
plot(t,p,'-b',[t(1) t(end)],[p_rms p_rms],'--r',[t(1) t(end)],[-p_rms -p_rms],'--r','Linewidth',2)
grid on
set(gca,'FontSize',20)
set(gcf, 'Renderer', 'zbuffer');
% ylim([-2000 2000])
xlim([t(1) t(end)])
xlabel('Время, с','FontSize',20)
ylabel('Избыточное давление, Па','FontSize',20)


%%
clc
close all
% -------------------------------------------------------------------------
% Фурье преобразование
phi_fft = fft(p,L)/L;
% забираем тольк N частот + постоянный член
y_fft = phi_fft(1:N+1);
% диапазон частот + постоянный член в 0
f_fft = (0:N)'*f_min_data;

% -----------------------------------
% in OpenFoam Prms_f = abs(y_fft); no need muliplacation by 2
% амплитуды (Fluent) равны размаху синуса , т.е. 2a. первый члена ряда a_0 (постоянный)
% не надо умножать на 2
% 
MagP_fft = 2*abs(y_fft);
MagP_fft(1) = abs(y_fft(1));

% мощность спектральной плотности  Power Spectral Density
PsdP_fft = 2*abs(y_fft).^2;
PsdP_fft(1) = abs(y_fft(1)).^2;
PsdP_fft = PsdP_fft*L*dt; % так как phi_fft уже поделено на L

% ------
% уровень звуковго давления -- Power Spectral Density in dB
% for Fluent it is SPL, in OPenFOAM is Power Spectral Density in dB
% https://www.sharcnet.ca/Software/Fluent6/html/ug/node1189.htm
p_ref=2e-5; % ссылочное давлени - порог слышимости на частоте 1 кГц
SPLp=10*log10(PsdP_fft/p_ref^2);


% -------------------------------------------------------------------------
% постоянный член (среднее значение давления) не отображаем
% так же ограничесваем отображение частот сверху и снизу
% f_min = f_min_data;
f_min = 20;
% f_max = 1/(2*dt);
f_max = 1e3;

% индексы
ind_f_min = find(f_fft>f_min,1);
ind_f_max = find(f_fft<f_max,1,'last');


% ----------------------------------------
% Сглаживающая функция
% SPLp_smooth = SPLp;

%  moving lowess loess rlowess rloess sgolay
span = 11;
% SPLp_smooth = smooth(SPLp,span,'sgolay',4);

% Apply 1/N-octave smoothing to a frequency spectrum
n_octave = 12;
PsdP_smooth = smoothSpectrum(PsdP_fft,f_fft,n_octave);
SPLp_smooth=10*log10(PsdP_smooth/p_ref^2);

% 1/n octave freq using base 2 filter frequencies
[fexact_l_c_u, fnominal_l_c_u, fnominal_str_l_c_u] = fract_oct_freq_band(n_octave, f_min, f_max, 0, 0,1000);

% convert narrowband data to 1/n octaveband data 
PsdP_N_octave = NarrowToNthOctave(PsdP_fft,f_fft,fexact_l_c_u);
indx = find(PsdP_N_octave<=0);
SPLp_N_octave=10*log10(PsdP_N_octave/p_ref^2);
SPLp_N_octave(indx) = 0;


% % Use Matlab Octave (use base 10 filter frequencies)
% [PsdP_N_octave_matlab,fexact_c_matlab] = poctave(PsdP_fft,fs,f_fft,'BandsPerOctave',n_octave,'FrequencyLimits',[f_min f_max],'psd');
% % [PsdP_N_octave_matlab,fexact_c_matlab] = poctave(p,fs,'BandsPerOctave',n_octave,'FrequencyLimits',[f_min f_max]);
% [fexact_l_c_u_malab, fnominal_l_c_u_malab, fnominal_str_l_c_u_malab] = fract_oct_freq_band(n_octave, f_min, f_max, 0, 1,1000);
% PsdP_N_octave_matlab = PsdP_N_octave_matlab./(fexact_l_c_u_malab(:,3) - fexact_l_c_u_malab(:,1));
% SPLp_N_octave_matlab=10*log10(PsdP_N_octave_matlab/p_ref^2);
% ----------------------

% semilogx(f_fft(ind_f_min:ind_f_max),MagP_fft(ind_f_min:ind_f_max),'-b','Linewidth',2)
% semilogx(f_fft(ind_f_min:ind_f_max),PsdP_fft((ind_f_min-1):(ind_f_max-1)),'-b+','Linewidth',2)
% semilogx(f_fft(ind_f_min:ind_f_max),SPLp(ind_f_min:ind_f_max),'-b','Linewidth',2)

hold on

semilogx(f_fft(ind_f_min:ind_f_max),SPLp(ind_f_min:ind_f_max),'-b','Linewidth',2)
semilogx(f_fft(ind_f_min:ind_f_max),SPLp_smooth(ind_f_min:ind_f_max),'-r','Linewidth',2)

hold off

semilogxOneNOctaveBar(SPLp_N_octave,fexact_l_c_u(:,1), fexact_l_c_u(:,3), 1,[1 0.5 1],0.5);
% semilogxOneNOctaveBar(SPLp_N_octave_matlab,fexact_l_c_u_malab(:,1), fexact_l_c_u(:,3), 1,[1 0.5 0.1],0.5);




% set(gca,'XScale','log')
% set(gca,'Xtick',-3:4); %// adjust manually; values in log scale
% set(gca,'Xticklabel',10.^get(gca,'Xtick'));

xlim([1e1 f_max])


% ------------------
% hold on 
% plot(f_fft(ind_f_min:ind_f_max),PsdP_fft((ind_f_min-1):(ind_f_max-1)),'-b','Linewidth',2)
% f_fft_1_3_start = one_third_freq_preferred(1:end-1) - (one_third_freq_preferred(2:end) - one_third_freq_preferred(1:end-1))/2; % не совсем верно, переделать
% stairs( f_fft_1_3_start(ind_f_min_1_3:ind_f_max_1_3),PsdP_1_3(ind_f_min_1_3:ind_f_max_1_3),'-r','Linewidth',2)
% plot(one_third_freq_preferred(ind_f_min_1_3:ind_f_max_1_3),PsdP_1_3(ind_f_min_1_3:ind_f_max_1_3),'-g+','Linewidth',2)
% hold off
% set(gca,'XScale','log')
% % ax = gca;
% % ax.XScale = 'log'; %make x-axis logarithmic
% % ax.XTick = one_third_freq_preferred; %set the tick values to be the center frequencies
% % ax.XTickLabelRotation = 90; %rotate the x-labels for better visibility
% % xlabel('Frequency in Hz');
% % ylabel('Sound pressure level in dB(A)');

legend('fft','smoosh specrum','1/12 -- octave band')
ylabel('Power Spectral Density in dB');

grid on
xlabel('Частота, Гц','FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Renderer', 'zbuffer');
% ylim([0 200])



return


%%
clc

% ------------------------------
% Постпроцессинг данных с Ensight
% data_file = 'x_0_1_y_-0_025_z_0_PSDf';
% data_file = 'x_0_15_y_0_z_0_025_PSDf';
data_file = 'e_24268_PSDf';
title_line = 1; % количество строк заголовка

fileID = fopen(data_file);

PSDf_ensight = textscan(fileID,repmat('%f',1,2), 'MultipleDelimsAsOne',true, 'HeaderLines',title_line);

f_fft_ensight(:,1) = PSDf_ensight{1};
PsdP_fft_ensight(:,1) = PSDf_ensight{2};

fclose(fileID);

% ------------------------------
SPLp_ensight=10*log10(PsdP_fft_ensight/p_ref^2);


% ------------------------------
hold on

semilogx(f_fft_ensight,SPLp_ensight,'-g','Linewidth',2)


hold off




% ------------------------------
xlim([1e1 f_max])
grid on
xlabel('Частота, Гц','FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Renderer', 'zbuffer');






%%
% Постпроцессинг вибро задачи
clc
close all
clear


addpath('mechanical_data')
% загрузка файла данных
data_file = 'p_x_0_5_y_0_025.txt';
% Количество строк в файле
fileID = fopen(data_file);
Cell_scan = textscan(fileID,'%s','Delimiter','\n');
num_line = length(Cell_scan{1,1});
title_line = 1; % количество строк заголовка
fclose(fileID);
clear Cell_scan
% считываем данные, без заголовка и последней строки 
u_over_f = dlmread(data_file,'\t',[title_line 0 (num_line - title_line) 1]);


semilogx(u_over_f(:,1),u_over_f(:,2),'-r+','Linewidth',2)
grid on
xlabel('Частота, Гц','FontSize',20)
ylabel('Амплитуда перемещений, м','FontSize',20)
set(gca,'FontSize',20)
set(gcf, 'Renderer', 'zbuffer');
ylim([0 3e-6])




