% clear all
close all
clc


%% TCC - SOLUÇÃO ANALÍTICA PARA MONOPOLO ESTACIONÁRIO:

%% DADOS UTILIZADOS NA SIMULAÇÃO

freq=100;         % frequência da onda 
S = 0.1;          % amplitude da vazão definida na simulação
t = 0.7;         % instante de tempo dos resultados
r_fonte = 0.05715/2; % raio da fonte

% CÁLCULO DA TEMPERATURA DA SIMULAÇÃO
c_0=340.29;       % velocidade do som utilizada no cálculo do comprimento de onda
c = 331.45;       % velocidade do som a 0 ºC
T0 = 273.15;      % 0 ºC em Kelvin
T = (c_0/c)^2*T0; % temperatura da simulação 

% CÁLCULO DA DENSIDADE
rho0 = 101325/(287.058*T); % equação dos gases ideais

% CÁLCULO DO TEMPO DE SIMULAÇÃO 
r_buffer = 204;                 % raio externo da zona de buffer
t_simulacao = (r_buffer)/c_0; % +6 para passar o primeiro pulso e fechar 210

%% SOLUÇÃO ANALÍTICA - AKHNOUKH

omega = freq*2*pi;
lambda = c_0/freq;
r = .0001:0.01:100;
r = .0001:0.01:130;

H_0 = besselh(0, (omega*r/c_0));       %Função de Hankel
G_t = (-1i*omega)*(1i/(4*c_0^2))*H_0;  %transformada de fourier da derivada dG/dt
P_2 = rho0*(c_0^2)*S*G_t;              %Amplitude considerando fonte pontual

p_2 = -1*imag(P_2*exp(1i*omega*t));    %calculado com fonte pontual

%% SOLUÇÃO ANALÍTICA - BUCKINGHAM

H_0_B = besselh(0, 2, (omega*r/c_0));  %Função de Hankel
P_2_B = -rho0*omega*S*H_0_B/4;
p_2_B = imag(P_2_B*exp(1i*omega*t));   %calculado com fonte pontual / solução invertida

%% SOLUÇÃO ANALÍTICA - JACOBSEN

Area = 2*pi*r_fonte;
Velocity = S/Area;
H_1_fonte_J = besselh(1, 2, (omega*r_fonte/c_0));

A0 = Velocity*1i*rho0*c_0/H_1_fonte_J;

H_0_J = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
P_2_J = A0*H_0_J;
p_2_J = imag(P_2_J*exp(1i*omega*t));                   %calculado com fonte pontual

%% RESULTADOS MONOPOLO ESTACIONÁRIO

Simulation8 = load('datawaveTransmissiveDamped_t070.txt');        % resultado da simulação
%Simulation16 = load('dataMalhaMenor_waveTransmissive_t040.txt');
% Simulation32 = load('datas01_ppw32_t050.txt');
% Simulation64 = load('datau01_ppw64_t050.txt');


figure(1)
hold on
%plot(r,p_2(1:length(r)), 'k-')                          % Analítico Akhnoukh com fonte pontual
%plot(r,p_2_B(1:length(r)), 'k-')                        % Analítico Buckingham com fonte pontual
plot(r/lambda,p_2_J(1:length(r)), 'k-')                        % Analítico Jacobsen com fonte pontual
plot(Simulation8(:,11)/lambda, (Simulation8(:,8))-101325, 'r-')   % Simulação
%plot(Simulation8(:,11)/lambda, (Simulation8(:,8))-101325, 'r-')   % damped
%plot(Simulation16(:,8)/lambda, (Simulation16(:,5))-101325, 'b--')   
%plot(Simulation16(:,11)/lambda, (Simulation16(:,8))-101325, 'b--')   % damped
%plot(Simulation32(:,8), (Simulation32(:,5))-101325)   
%plot(Simulation64(:,8), (Simulation64(:,5))-101325)   
%plot(Simulation(:,13), (Simulation(:,8))-101325, 'r-') % Simulação com
%meta data
hold off
grid on
xlabel('r/\lambda')
ylabel('p (Pa)')
legend('Analítico','waveTransmissive', 'Location','best');
set(gca,'FontSize',20)


%% ERRO REFINO

% %Simulation8 = load('datas01_8ppw_t050.txt');        % resultado da simulação
% Simulation16 = load('datas01_ppw16_t050.txt');
% Simulation32 = load('datas01_ppw32_t050_2.txt');
% Simulation64 = load('datas01_ppw64_t050.txt');
% 
% % %% 8 PPW
% % % Analítico
% % r_erro8 = Simulation8(:,8);
% % H_0_E8 = besselh(0, 2, (omega*r_erro8/c_0));                  %Função de Hankel
% % P_2_E8 = A0*H_0_E8;
% % p_2_E8 = imag(P_2_E8*exp(1i*omega*t));                   %calculado com fonte pontual
% % 
% % % RMS
% % Num1_8 = (Simulation8(:,5)-101325)-p_2_E8;
% % Num2_8 = Num1_8.^2;
% % Num_8 = trapz(r_erro8,Num2_8);
% % Den_8 = trapz(r_erro8,p_2_E8.^2);
% % RMS_8 = Num_8/Den_8;
% 
% %% 16 PPW
% % Analítico
% r_erro16 = Simulation16(:,8);
% H_0_E16 = besselh(0, 2, (omega*r_erro16/c_0));                  %Função de Hankel
% P_2_E16 = A0*H_0_E16;
% p_2_E16 = imag(P_2_E16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_16 = (Simulation16(:,5)-101325)-p_2_E16;
% Num2_16 = Num1_16.^2;
% Num_16 = trapz(r_erro16,Num2_16);
% Den_16 = trapz(r_erro16,p_2_E16.^2);
% RMS_16 = Num_16/Den_16;
% 
% %% 32 PPW
% % Analítico
% r_erro32 = Simulation32(:,8);
% H_0_E32 = besselh(0, 2, (omega*r_erro32/c_0));                  %Função de Hankel
% P_2_E32 = A0*H_0_E32;
% p_2_E32 = imag(P_2_E32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_32 = (Simulation32(:,5)-101325)-p_2_E32;
% Num2_32 = Num1_32.^2;
% Num_32 = trapz(r_erro32,Num2_32);
% Den_32 = trapz(r_erro32,p_2_E32.^2);
% RMS_32 = Num_32/Den_32;
% 
% %% 64 PPW
% % Analítico
% r_erro64 = Simulation64(:,8);
% H_0_E64 = besselh(0, 2, (omega*r_erro64/c_0));                  %Função de Hankel
% P_2_E64 = A0*H_0_E64;
% p_2_E64 = imag(P_2_E64*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_64 = (Simulation64(:,5)-101325)-p_2_E64;
% Num2_64 = Num1_64.^2;
% Num_64 = trapz(r_erro64,Num2_64);
% Den_64 = trapz(r_erro64,p_2_E64.^2);
% RMS_64 = Num_64/Den_64;
% 
% %% PPW X RMS
% 
% RMS = [RMS_16, RMS_32, RMS_64]*100;
% PPW = [16, 32, 64];
% 
% figure(2)
% hold on
% plot(PPW,RMS,'-o','color',[0.6350 0.0780 0.1840])
% hold off
% xticks([16 32 64]);
% xlabel('PPW')
% ylabel('Erro (%)')
% set(gca,'FontSize',15)
% % for i = 1:length(PPW)
% %     text(PPW(i)+0.5,RMS(i)+0.5,['(',num2str(PPW(i)),',',num2str(RMS(i)),')'])
% % end

%% ZONA ABSORÇÃO

% Simulation8 = load('datawaveTransmissiveDamped_t070.txt');        % resultado da simulação
%  Simulation16 = load('datafixedValue_damped_t070.txt');
%  Simulation32 = load('dataMalhaMenor_waveTransmissiveDamped_t040.txt');
%   Simulation64 = load('dataMalhaMenor_fixedValueDamped_t040.txt');
% 
% 
% %% WAVETRANSMISSIVE MALHA 3
% % Analítico
% r_erro8 = Simulation8(:,11);
% H_0_E8 = besselh(0, 2, (omega*r_erro8/c_0));                  %Função de Hankel
% P_2_E8 = A0*H_0_E8;
% p_2_E8 = imag(P_2_E8*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_8 = (Simulation8(:,8)-101325)-p_2_E8;
% Num2_8 = Num1_8.^2;
% Num_8 = trapz(r_erro8,Num2_8);
% Den_8 = trapz(r_erro8,p_2_E8.^2);
% RMS_8 = Num_8/Den_8;
% 
% %% FIXED VALUE MALHA 3
% % Analítico
% r_erro16 = Simulation16(:,11);
% H_0_E16 = besselh(0, 2, (omega*r_erro16/c_0));                  %Função de Hankel
% P_2_E16 = A0*H_0_E16;
% p_2_E16 = imag(P_2_E16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_16 = (Simulation16(:,8)-101325)-p_2_E16;
% Num2_16 = Num1_16.^2;
% Num_16 = trapz(r_erro16,Num2_16);
% Den_16 = trapz(r_erro16,p_2_E16.^2);
% RMS_16 = Num_16/Den_16;
% 
% %% WAVE TRANSMISSIVE MALHA 4
% t = 0.4;
% % Analítico
% r_erro32 = Simulation32(:,11);
% H_0_E32 = besselh(0, 2, (omega*r_erro32/c_0));                  %Função de Hankel
% P_2_E32 = A0*H_0_E32;
% p_2_E32 = imag(P_2_E32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_32 = (Simulation32(:,8)-101325)-p_2_E32;
% Num2_32 = Num1_32.^2;
% Num_32 = trapz(r_erro32,Num2_32);
% Den_32 = trapz(r_erro32,p_2_E32.^2);
% RMS_32 = Num_32/Den_32;
% 
% %% FIXED VALUE MALHA 4
% % Analítico
% r_erro64 = Simulation64(:,11);
% H_0_E64 = besselh(0, 2, (omega*r_erro64/c_0));                  %Função de Hankel
% P_2_E64 = A0*H_0_E64;
% p_2_E64 = imag(P_2_E64*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_64 = (Simulation64(:,8)-101325)-p_2_E64;
% Num2_64 = Num1_64.^2;
% Num_64 = trapz(r_erro64,Num2_64);
% Den_64 = trapz(r_erro64,p_2_E64.^2);
% RMS_64 = Num_64/Den_64;
% 
% 
% %% PPW X RMS
% 
% RMS = [RMS_8, RMS_16, RMS_32, RMS_64]*100
% PPW = [8, 16, 32, 64];
% 
% % figure(2)
% % hold on
% % plot(PPW,RMS,'-o','color',[0.6350 0.0780 0.1840])
% % hold off
% % xticks([8 16 32 64]);
% % xlabel('PPW')
% % ylabel('Erro (%)')
% % % for i = 1:length(PPW)
% % %     text(PPW(i)+0.5,RMS(i)+0.5,['(',num2str(PPW(i)),',',num2str(RMS(i)),')'])
% % % end

%% ERRO CONDIÇÃO DE CONTORNO

% Simulation8 = load('dataMalhaMenor_waveTransmissive_t040.txt');        % resultado da simulação
%  Simulation16 = load('datas01_ppw32_t040.txt');
%  Simulation32 = load('dataMalhaMenor_fixedValue_t040.txt');
% 
% %% WAVETRANSMISSIVE
% % Analítico
% r_erro8 = Simulation8(:,8);
% H_0_E8 = besselh(0, 2, (omega*r_erro8/c_0));                  %Função de Hankel
% P_2_E8 = A0*H_0_E8;
% p_2_E8 = imag(P_2_E8*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_8 = (Simulation8(:,5)-101325)-p_2_E8;
% Num2_8 = Num1_8.^2;
% Num_8 = trapz(r_erro8,Num2_8);
% Den_8 = trapz(r_erro8,p_2_E8.^2);
% RMS_8 = Num_8/Den_8;
% 
% %% ZONA SAÍDA
% % Analítico
% r_erro16 = Simulation16(:,8);
% H_0_E16 = besselh(0, 2, (omega*r_erro16/c_0));                  %Função de Hankel
% P_2_E16 = A0*H_0_E16;
% p_2_E16 = imag(P_2_E16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_16 = (Simulation16(:,5)-101325)-p_2_E16;
% Num2_16 = Num1_16.^2;
% Num_16 = trapz(r_erro16,Num2_16);
% Den_16 = trapz(r_erro16,p_2_E16.^2);
% RMS_16 = Num_16/Den_16;
% 
% %% FIXED VALUE
% % Analítico
% r_erro32 = Simulation32(:,8);
% H_0_E32 = besselh(0, 2, (omega*r_erro32/c_0));                  %Função de Hankel
% P_2_E32 = A0*H_0_E32;
% p_2_E32 = imag(P_2_E32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_32 = (Simulation32(:,5)-101325)-p_2_E32;
% Num2_32 = Num1_32.^2;
% Num_32 = trapz(r_erro32,Num2_32);
% Den_32 = trapz(r_erro32,p_2_E32.^2);
% RMS_32 = Num_32/Den_32;
% 
% 
% %% PPW X RMS
% 
% RMS = [RMS_8, RMS_16, RMS_32]*100
% PPW = [8, 16, 32];
% 
% % figure(2)
% % hold on
% % plot(PPW,RMS,'-o','color',[0.6350 0.0780 0.1840])
% % hold off
% % xticks([8 16 32]);
% % xlabel('PPW')
% % ylabel('Erro (%)')
% % % for i = 1:length(PPW)
% % %     text(PPW(i)+0.5,RMS(i)+0.5,['(',num2str(PPW(i)),',',num2str(RMS(i)),')'])
% % % end

%% ERRO TEMPORAL

% Euler10 = load('dataEuler10_t050.txt');        % resultado da simulação
% Euler20 = load('dataEuler20_t050.txt');
% Euler40 = load('dataEuler40_t050.txt');
% Euler80 = load('dataEuler80_t050.txt');
% Euler4000 = load('dataEuler80_t050.txt');
% 
% Backward10 = load('dataBackward10_t050.txt');        % resultado da simulação
% Backward20 = load('dataBackward20_t050.txt');
% Backward40 = load('dataBackward40_t050.txt');
% Backward80 = load('dataBackward80_t050.txt');
% Backward4000 = load('datas01_ppw32_t050_2.txt');
% 
% %% 10 DIVISOES
% % Analítico
% r_erroEuler10 = Euler10(:,8);
% H_0_EEuler10 = besselh(0, 2, (omega*r_erroEuler10/c_0));                  %Função de Hankel
% P_2_EEuler10 = A0*H_0_EEuler10;
% p_2_EEuler10 = imag(P_2_EEuler10*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroBackward10 = Backward10(:,8);
% H_0_EBackward10 = besselh(0, 2, (omega*r_erroBackward10/c_0));                  %Função de Hankel
% P_2_EBackward10 = A0*H_0_EBackward10;
% p_2_EBackward10 = imag(P_2_EBackward10*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Euler10 = (Euler10(:,5)-101325)-p_2_EEuler10;
% Num2_Euler10 = Num1_Euler10.^2;
% Num_Euler10 = trapz(r_erroEuler10,Num2_Euler10);
% Den_Euler10 = trapz(r_erroEuler10,p_2_EEuler10.^2);
% RMS_Euler10 = Num_Euler10/Den_Euler10;
% 
% Num1_Backward10 = (Backward10(:,5)-101325)-p_2_EBackward10;
% Num2_Backward10 = Num1_Backward10.^2;
% Num_Backward10 = trapz(r_erroBackward10,Num2_Backward10);
% Den_Backward10 = trapz(r_erroBackward10,p_2_EBackward10.^2);
% RMS_Backward10 = Num_Backward10/Den_Backward10;
% 
% %% 20 DIVISOES
% % Analítico
% r_erroEuler20 = Euler20(:,8);
% H_0_EEuler20 = besselh(0, 2, (omega*r_erroEuler20/c_0));                  %Função de Hankel
% P_2_EEuler20 = A0*H_0_EEuler20;
% p_2_EEuler20 = imag(P_2_EEuler20*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroBackward20 = Backward20(:,8);
% H_0_EBackward20 = besselh(0, 2, (omega*r_erroBackward20/c_0));                  %Função de Hankel
% P_2_EBackward20 = A0*H_0_EBackward20;
% p_2_EBackward20 = imag(P_2_EBackward20*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Euler20 = (Euler20(:,5)-101325)-p_2_EEuler20;
% Num2_Euler20 = Num1_Euler20.^2;
% Num_Euler20 = trapz(r_erroEuler20,Num2_Euler20);
% Den_Euler20 = trapz(r_erroEuler20,p_2_EEuler20.^2);
% RMS_Euler20 = Num_Euler20/Den_Euler20;
% 
% Num1_Backward20 = (Backward20(:,5)-101325)-p_2_EBackward20;
% Num2_Backward20 = Num1_Backward20.^2;
% Num_Backward20 = trapz(r_erroBackward20,Num2_Backward20);
% Den_Backward20 = trapz(r_erroBackward20,p_2_EBackward20.^2);
% RMS_Backward20 = Num_Backward20/Den_Backward20;
% 
% %% 40 DIVISOES
% % Analítico
% r_erroEuler40 = Euler40(:,8);
% H_0_EEuler40 = besselh(0, 2, (omega*r_erroEuler40/c_0));                  %Função de Hankel
% P_2_EEuler40 = A0*H_0_EEuler40;
% p_2_EEuler40 = imag(P_2_EEuler40*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroBackward40 = Backward40(:,8);
% H_0_EBackward40 = besselh(0, 2, (omega*r_erroBackward40/c_0));                  %Função de Hankel
% P_2_EBackward40 = A0*H_0_EBackward40;
% p_2_EBackward40 = imag(P_2_EBackward40*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Euler40 = (Euler40(:,5)-101325)-p_2_EEuler40;
% Num2_Euler40 = Num1_Euler40.^2;
% Num_Euler40 = trapz(r_erroEuler40,Num2_Euler40);
% Den_Euler40 = trapz(r_erroEuler40,p_2_EEuler40.^2);
% RMS_Euler40 = Num_Euler40/Den_Euler40;
% 
% Num1_Backward40 = (Backward40(:,5)-101325)-p_2_EBackward40;
% Num2_Backward40 = Num1_Backward40.^2;
% Num_Backward40 = trapz(r_erroBackward40,Num2_Backward40);
% Den_Backward40 = trapz(r_erroBackward40,p_2_EBackward40.^2);
% RMS_Backward40 = Num_Backward40/Den_Backward40;
% 
% %% 80 DIVISOES
% % Analítico
% r_erroEuler80 = Euler80(:,8);
% H_0_EEuler80 = besselh(0, 2, (omega*r_erroEuler80/c_0));                  %Função de Hankel
% P_2_EEuler80 = A0*H_0_EEuler80;
% p_2_EEuler80 = imag(P_2_EEuler80*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroBackward80 = Backward80(:,8);
% H_0_EBackward80 = besselh(0, 2, (omega*r_erroBackward80/c_0));                  %Função de Hankel
% P_2_EBackward80 = A0*H_0_EBackward80;
% p_2_EBackward80 = imag(P_2_EBackward80*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Euler80 = (Euler80(:,5)-101325)-p_2_EEuler80;
% Num2_Euler80 = Num1_Euler80.^2;
% Num_Euler80 = trapz(r_erroEuler80,Num2_Euler80);
% Den_Euler80 = trapz(r_erroEuler80,p_2_EEuler80.^2);
% RMS_Euler80 = Num_Euler80/Den_Euler80;
% 
% Num1_Backward80 = (Backward80(:,5)-101325)-p_2_EBackward80;
% Num2_Backward80 = Num1_Backward80.^2;
% Num_Backward80 = trapz(r_erroBackward80,Num2_Backward80);
% Den_Backward80 = trapz(r_erroBackward80,p_2_EBackward80.^2);
% RMS_Backward80 = Num_Backward80/Den_Backward80;
% 
% %% 4000 DIVISOES
% % Analítico
% r_erroEuler4000 = Euler4000(:,8);
% H_0_EEuler4000 = besselh(0, 2, (omega*r_erroEuler4000/c_0));                  %Função de Hankel
% P_2_EEuler4000 = A0*H_0_EEuler4000;
% p_2_EEuler4000 = imag(P_2_EEuler4000*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroBackward4000 = Backward4000(:,8);
% H_0_EBackward4000 = besselh(0, 2, (omega*r_erroBackward4000/c_0));                  %Função de Hankel
% P_2_EBackward4000 = A0*H_0_EBackward4000;
% p_2_EBackward4000 = imag(P_2_EBackward4000*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Euler4000 = (Euler4000(:,5)-101325)-p_2_EEuler4000;
% Num2_Euler4000 = Num1_Euler4000.^2;
% Num_Euler4000 = trapz(r_erroEuler4000,Num2_Euler4000);
% Den_Euler4000 = trapz(r_erroEuler4000,p_2_EEuler4000.^2);
% RMS_Euler4000 = Num_Euler4000/Den_Euler4000;
% 
% Num1_Backward4000 = (Backward4000(:,5)-101325)-p_2_EBackward4000;
% Num2_Backward4000 = Num1_Backward4000.^2;
% Num_Backward4000 = trapz(r_erroBackward4000,Num2_Backward4000);
% Den_Backward4000 = trapz(r_erroBackward4000,p_2_EBackward4000.^2);
% RMS_Backward4000 = Num_Backward4000/Den_Backward4000;
% 
% %% n X RMS
% 
% RMS_Euler = [RMS_Euler10, RMS_Euler20, RMS_Euler40, RMS_Euler80]*100;
% RMS_Backward = [RMS_Backward10, RMS_Backward20, RMS_Backward40, RMS_Backward80]*100;
% n = [10, 20, 40, 80];
% 
% figure(3)
% hold on
% plot(n,RMS_Euler,'-o','color','k')
% plot(n,RMS_Backward,'--o','color','b')
% hold off
% xticks([10 20 40 80]);
% xlabel('n')
% ylabel('Erro (%)')
% legend('Euler', 'Backward', 'Location','best');

%% ERRO ESQUEMAS

% Upwind16 = load('dataUpwind_16ppw_t050.txt');        % resultado da simulação
% linearUpwind16 = load('dataLinearUpwind_16ppw_t050.txt');
% limitedLinear16 = load('dataLimited_16ppw_t050.txt');
% vanLeer16 = load('datas01_ppw16_t050.txt');
% 
% Upwind32 = load('dataUpwind_32ppw_t050.txt');        % resultado da simulação
% linearUpwind32 = load('dataLinearUpwind_32ppw_t050.txt');
% limitedLinear32 = load('dataLimited_32ppw_t050.txt');
% vanLeer32 = load('datas01_ppw32_t050.txt');
% 
% %% UPWIND
% % Analítico
% r_erroUpwind16 = Upwind16(:,8);
% H_0_EUpwind16 = besselh(0, 2, (omega*r_erroUpwind16/c_0));                  %Função de Hankel
% P_2_EUpwind16 = A0*H_0_EUpwind16;
% p_2_EUpwind16 = imag(P_2_EUpwind16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_erroUpwind32 = Upwind32(:,8);
% H_0_EUpwind32 = besselh(0, 2, (omega*r_erroUpwind32/c_0));                  %Função de Hankel
% P_2_EUpwind32 = A0*H_0_EUpwind32;
% p_2_EUpwind32 = imag(P_2_EUpwind32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_Upwind16 = (Upwind16(:,5)-101325)-p_2_EUpwind16;
% Num2_Upwind16 = Num1_Upwind16.^2;
% Num_Upwind16 = trapz(r_erroUpwind16,Num2_Upwind16);
% Den_Upwind16 = trapz(r_erroUpwind16,p_2_EUpwind16.^2);
% RMS_Upwind16 = Num_Upwind16/Den_Upwind16;
% 
% Num1_Upwind32 = (Upwind32(:,5)-101325)-p_2_EUpwind32;
% Num2_Upwind32 = Num1_Upwind32.^2;
% Num_Upwind32 = trapz(r_erroUpwind32,Num2_Upwind32);
% Den_Upwind32 = trapz(r_erroUpwind32,p_2_EUpwind32.^2);
% RMS_Upwind32 = Num_Upwind32/Den_Upwind32;
% 
% %% LINEAR UPWIND
% % Analítico
% r_errolinearUpwind16 = linearUpwind16(:,8);
% H_0_ElinearUpwind16 = besselh(0, 2, (omega*r_errolinearUpwind16/c_0));                  %Função de Hankel
% P_2_ElinearUpwind16 = A0*H_0_ElinearUpwind16;
% p_2_ElinearUpwind16 = imag(P_2_ElinearUpwind16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_errolinearUpwind32 = linearUpwind32(:,8);
% H_0_ElinearUpwind32 = besselh(0, 2, (omega*r_errolinearUpwind32/c_0));                  %Função de Hankel
% P_2_ElinearUpwind32 = A0*H_0_ElinearUpwind32;
% p_2_ElinearUpwind32 = imag(P_2_ElinearUpwind32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_linearUpwind16 = (linearUpwind16(:,5)-101325)-p_2_ElinearUpwind16;
% Num2_linearUpwind16 = Num1_linearUpwind16.^2;
% Num_linearUpwind16 = trapz(r_errolinearUpwind16,Num2_linearUpwind16);
% Den_linearUpwind16 = trapz(r_errolinearUpwind16,p_2_ElinearUpwind16.^2);
% RMS_linearUpwind16 = Num_linearUpwind16/Den_linearUpwind16;
% 
% Num1_linearUpwind32 = (linearUpwind32(:,5)-101325)-p_2_ElinearUpwind32;
% Num2_linearUpwind32 = Num1_linearUpwind32.^2;
% Num_linearUpwind32 = trapz(r_errolinearUpwind32,Num2_linearUpwind32);
% Den_linearUpwind32 = trapz(r_errolinearUpwind32,p_2_ElinearUpwind32.^2);
% RMS_linearUpwind32 = Num_linearUpwind32/Den_linearUpwind32;
% 
% %% LIMITED LINEAR
% % Analítico
% r_errolimitedLinear16 = limitedLinear16(:,8);
% H_0_ElimitedLinear16 = besselh(0, 2, (omega*r_errolimitedLinear16/c_0));                  %Função de Hankel
% P_2_ElimitedLinear16 = A0*H_0_ElimitedLinear16;
% p_2_ElimitedLinear16 = imag(P_2_ElimitedLinear16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_errolimitedLinear32 = limitedLinear32(:,8);
% H_0_ElimitedLinear32 = besselh(0, 2, (omega*r_errolimitedLinear32/c_0));                  %Função de Hankel
% P_2_ElimitedLinear32 = A0*H_0_ElimitedLinear32;
% p_2_ElimitedLinear32 = imag(P_2_ElimitedLinear32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_limitedLinear16 = (limitedLinear16(:,5)-101325)-p_2_ElimitedLinear16;
% Num2_limitedLinear16 = Num1_limitedLinear16.^2;
% Num_limitedLinear16 = trapz(r_errolimitedLinear16,Num2_limitedLinear16);
% Den_limitedLinear16 = trapz(r_errolimitedLinear16,p_2_ElimitedLinear16.^2);
% RMS_limitedLinear16 = Num_limitedLinear16/Den_limitedLinear16;
% 
% Num1_limitedLinear32 = (limitedLinear32(:,5)-101325)-p_2_ElimitedLinear32;
% Num2_limitedLinear32 = Num1_limitedLinear32.^2;
% Num_limitedLinear32 = trapz(r_errolimitedLinear32,Num2_limitedLinear32);
% Den_limitedLinear32 = trapz(r_errolimitedLinear32,p_2_ElimitedLinear32.^2);
% RMS_limitedLinear32 = Num_limitedLinear32/Den_limitedLinear32;
% 
% %% VAN LEER
% % Analítico
% r_errovanLeer16 = vanLeer16(:,8);
% H_0_EvanLeer16 = besselh(0, 2, (omega*r_errovanLeer16/c_0));                  %Função de Hankel
% P_2_EvanLeer16 = A0*H_0_EvanLeer16;
% p_2_EvanLeer16 = imag(P_2_EvanLeer16*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% r_errovanLeer32 = vanLeer32(:,8);
% H_0_EvanLeer32 = besselh(0, 2, (omega*r_errovanLeer32/c_0));                  %Função de Hankel
% P_2_EvanLeer32 = A0*H_0_EvanLeer32;
% p_2_EvanLeer32 = imag(P_2_EvanLeer32*exp(1i*omega*t));                   %calculado com fonte pontual
% 
% % RMS
% Num1_vanLeer16 = (vanLeer16(:,5)-101325)-p_2_EvanLeer16;
% Num2_vanLeer16 = Num1_vanLeer16.^2;
% Num_vanLeer16 = trapz(r_errovanLeer16,Num2_vanLeer16);
% Den_vanLeer16 = trapz(r_errovanLeer16,p_2_EvanLeer16.^2);
% RMS_vanLeer16 = Num_vanLeer16/Den_vanLeer16;
% 
% Num1_vanLeer32 = (vanLeer32(:,5)-101325)-p_2_EvanLeer32;
% Num2_vanLeer32 = Num1_vanLeer32.^2;
% Num_vanLeer32 = trapz(r_errovanLeer32,Num2_vanLeer32);
% Den_vanLeer32 = trapz(r_errovanLeer32,p_2_EvanLeer32.^2);
% RMS_vanLeer32 = Num_vanLeer32/Den_vanLeer32;
% 
% %% ESQUEMA X RMS
% 
% RMS_esquema16 = [RMS_limitedLinear16, RMS_linearUpwind16, RMS_Upwind16, RMS_vanLeer16]*100;
% RMS_esquema32 = [RMS_limitedLinear32, RMS_linearUpwind32, RMS_Upwind32, RMS_vanLeer32]*100;
% 
% RMS_limitedLinear = [RMS_limitedLinear16, RMS_limitedLinear32]*100;
% RMS_linearUpwind = [RMS_linearUpwind16, RMS_linearUpwind32]*100;
% RMS_Upwind = [RMS_Upwind16, RMS_Upwind32]*100;
% RMS_vanLeer = [RMS_vanLeer16, RMS_vanLeer32]*100;
% 
% esquemas_PPW = [16, 32];
% 
% esquemas_label16=["limitedLinear","linearUpwind","upwind","vanLeer"];
% esquemas16 = categorical(esquemas_label16);
% 
% % figure(3)
% % hold on
% % plot(esquemas16,RMS_esquema16,'-o','color','b')
% % hold off
% % xlabel('Esquemas Numericos')
% % ylabel('RMS (%)')
% % legend('16 PPW', 'Location','best');
% % 
% esquemas_label32=["limitedLinear","linearUpwind","upwind","vanLeer"];
% esquemas32 = categorical(esquemas_label32);
% % 
% % figure(4)
% % hold on
% % plot(esquemas32,RMS_esquema32,'-o','color','b')
% % hold off
% % xlabel('Esquemas Numericos')
% % ylabel('Erro (%)')
% % legend('32 PPW', 'Location','best');
% 
% figure(5)
% hold on
% bar(esquemas16,RMS_esquema16,'c')
% bar(esquemas32,RMS_esquema32,'r')
% hold off
% grid on
% xlabel('Esquemas Numericos')
% ylabel('Erro (%)')
% legend('16 PPW', '32 PPW', 'Location','best');
% set(gca,'FontSize',20)
% 
% % figure(6)
% % hold on
% % plot(esquemas_PPW,RMS_limitedLinear,'-o','color','r')
% % plot(esquemas_PPW,RMS_linearUpwind,'-*','color','k')
% % plot(esquemas_PPW,RMS_Upwind,'-x','color','b')
% % plot(esquemas_PPW,RMS_vanLeer,'-s','color','m')
% % hold off
% % %xticks([16 32]);
% % xlabel('PPW')
% % ylabel('Erro (%)')
% % legend('limitedLinear','linearUpwind', 'upwind', 'vanLeer', 'Location','best');
% % set(gca,'FontSize',20)


%% MONOPOLO COM ESCOAMENTO


% %% SOLUÇÃO ANALÍTICA
% 
% M = 0.1;
% U = M*c_0;
% y = 0;
% x = -100:0.13:100;
% k = omega/c_0; % número de onda
% Gx = zeros(1,length(x)); 
% Gt_escoamento = zeros(1,length(x));
% 
% ksi = omega*sqrt(x.^2+(1-M^2)*y^2)/((1-M^2)*c_0);
% eta = -1i*M/(1-M^2)*k*x-1i*omega*t;
% H_0_escoamento = besselh(0, ksi);
% H_1_escoamento = besselh(1, ksi);
% 
% % DERIVADA DA FUNÇÃO GREEN EM X
% for i=1:length(x)
% Gx(i) = omega/(4*c_0^3*(1-M^2)^(3/2))*(M*H_0_escoamento(i)-1i*x(i)*H_1_escoamento(i)/sqrt(x(i)^2+(1-M^2)*y^2))*exp(eta(i));
% end
% 
% % DERIVADA DA FUNÇÃO GREEN EM t
% for i=1:length(x)
% Gt_escoamento(i) = 1i/(4*c_0^2*sqrt(1-M^2))*H_0_escoamento(i)*exp(eta(i))*(-1i*omega);
% end
% 
% p_flow = -imag(rho0*(c_0^2)*S*(Gt_escoamento+M*Gx));
% 
% %% RESULTADOS - MONOPOLO COM ESCOAMENTO
% 
% Escoamento = load('dataM01_s10_32ppwref_t0343.txt');        % resultado da simulação
% 
% figure(2)
% hold on
% plot(x,p_flow, 'k-')                                  % Analítico Akhnoukh com escoamento
% plot(Escoamento(:,8)-100, (Escoamento(:,5))-101325, 'r--')   % Simulação
% %plot(Escoamento(:,8), (Escoamento(:,5))-101325, 'r-')   % Simulação
% %plot(Escoamento(:,11), (Escoamento(:,8))-101325, 'r-')   % damped
% hold off
% xlabel('x (m)')
% ylabel('p (Pa)')
% legend('Analítico','Numérico', 'Location','best');
% 
% %sqrt((pi*(204^2-r_fonte^2))/394309)