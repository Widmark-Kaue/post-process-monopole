%clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TCC - SOLUÇÃO ANALÍTICA PARA MONOPOLO ESTACIONÁRIO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DADOS UTILIZADOS NA SIMULAÇÃO

%c_0=1;           % velocidade do som para fechar lambda=30
freq=100;         % frequência da onda 
P_sim=10;         % amplitude da pressão definida na simulação
t = 0.17;         % instante de tempo dos resultados

% CÁLCULO DA TEMPERATURA DA SIMULAÇÃO
c_0=340.29;       % velocidade do som utilizada no cálculo do comprimento de onda
c = 331.45;       % velocidade do som a 0 ºC
T0 = 273.15;      % 0 ºC em Kelvin
T = (c_0/c)^2*T0; % temperatura da simulação 
%omega=pi()/15;   % Dissertação
%A=23.206;        % Constante para fonte pontual - definir a partir da aplicação de condição de contorno
%A=1;

% CÁLCULO DA DENSIDADE
rho0 = 101325/(287.058*T); % equação dos gases ideais

% CÁLCULO DO TEMPO DE SIMULAÇÃO 
r_buffer = 204;                 % raio externo da zona de buffer
t_simulacao = (r_buffer+6)/c_0; % +6 para passar o primeiro pulso e fechar 210

%% SOLUÇÃO ANALÍTICA
%freq=omega/(2*pi());
omega=freq*2*pi;
lambda=c_0/freq;
r=.0001:0.1:100;

%solução analitica da dissertação
epsilon=1;
alpha=log(2)/(9);
f=epsilon*exp(-alpha*r.^2);


H_0=besselh(0, (omega*r/c_0));      %Função de Hankel
G_t=(-1i*omega)*(1i/(4*c_0^2))*H_0; %transformada de fourier da derivada dG/dt

%convolução via FFT espacial
%f_k=fft(f);
%G_k=fft(G_t);
%P=ifft(f_k.*G_k);

P=conv(f,G_t); %convolução

%% CÁLCULO DA CONSTANTE "A" PELA CONDIÇÃO DE CONTORNO

r_fonte = 0.05715/2;
H_0_fonte=besselh(0, (omega*r_fonte/c_0));      %Função de Hankel
G_t_fonte=(-1i*omega)*(1i/(4*c_0^2))*H_0_fonte; %transformada de fourier da derivada dG/dt

A = abs(P_sim/G_t_fonte);
P_2=A*G_t;                                     %Amplitude considerando fonte pontual

p_2=-1*imag(P_2*exp(-1i*omega*t));              %calculado com fonte pontual

%% SOLUÇÃO ANALÍTICA - BUCKINGHAM

H_0_B=besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
% UTILIZANDO PARÂMETRO DE PRESSÃO (COMENTAR A SOLUÇÃO DA VELOCIDADE PARA
% CALCULAR)

H_0_fonte_B=besselh(0, 2, (omega*r_fonte/c_0));      %Função de Hankel
A_B = abs(P_sim/H_0_fonte_B);
P_2_B=A_B*H_0_B;
p_2_B = imag(P_2_B*exp(1i*omega*t)); %calculado com fonte pontual

% UTILIZANDO PARÂMETRO DE VELOCIDADE
S = 0.1;                              % Vazão volumétrica definida na simulação por unidade de espessura do domínio
P_2_B = -rho0*omega*S*H_0_B/4;
p_2_B = imag(P_2_B*exp(1i*omega*t)); %calculado com fonte pontual / solução invertida

%% SOLUÇÃO ANALÍTICA - JACOBSEN

Area = 2*pi*r_fonte;
Velocity = S/Area;
H_1_fonte_J=besselh(1, 2, (omega*r_fonte/c_0));

A0 = Velocity*1i*rho0*c_0/H_1_fonte_J;

H_0_J=besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
P_2_J = A0*H_0_J;
p_2_J = imag(P_2_J*exp(1i*omega*t)); %calculado com fonte pontual

%% RESULTADOS MONOPOLO ESTACIONÁRIO
%p_=-1*imag(P*exp(-1i*omega*t));         % Aqui utilizei -Im pois a excitação é um seno na dissertação

% Thesis=load('akhnoukh.txt');           % analítico da dissertação extraído do gráfico
Simulation = load('datau01_ppw32_t017.txt');        % resultado da simulação

% figure(1)
% %plot(r,p_(1:length(r)))                                % Plot do analítico com fonte da dissertação
% % plot(Thesis(:,1), Thesis(:,2), 'r-')                  % Analítico extraído do gráfico
% hold on
% %plot(r,p_2(1:length(r)), 'k-')                          % Analítico Akhnoukh com fonte pontual
% %plot(r,p_2_B(1:length(r)), 'k-')                        % Analítico Buckingham com fonte pontual
% plot(r,p_2_J(1:length(r)), 'k-')                        % Analítico Jacobsen com fonte pontual
% plot(Simulation(:,8), (Simulation(:,5))-101325, 'r-')   % Simulação
% %plot(Simulation(:,13), (Simulation(:,8))-101325, 'r-') % Simulação com
% %meta data
% hold off
% xlabel('r (m)')
% ylabel('p (Pa)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MONOPOLO COM ESCOAMENTO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOLUÇÃO ANALÍTICA
M = 0.3;
y = 0;
x = -100:0.13:100;
k = omega/c_0; % número de onda
Gx = zeros(1,length(x)); 
Gt_escoamento = zeros(1,length(x));

ksi = omega*sqrt(x.^2+(1-M^2)*y^2)/((1-M^2)*c_0);
eta = -1i*M/(1-M^2)*k*x-1i*omega*t;
H_0_escoamento = besselh(0, ksi);
H_1_escoamento = besselh(1, ksi);

% DERIVADA DA FUNÇÃO GREEN EM X
for i=1:length(x)
Gx(i) = omega/(4*c_0^3*(1-M^2)^(3/2))*(M*H_0_escoamento(i)-1i*x(i)*H_1_escoamento(i)/sqrt(x(i)^2+(1-M^2)*y^2))*exp(eta(i));
end

% DERIVADA DA FUNÇÃO GREEN EM t
for i=1:length(x)
Gt_escoamento(i) = eta(i)*1i/(4*c_0^2*sqrt(1-M^2))*H_0_escoamento(i)*exp(eta(i))*(-1i*omega);
end

p_flow=-imag(rho0*(c_0^2)*S*(Gt_escoamento+M*Gx));

%% RESULTADOS - MONOPOLO COM ESCOAMENTO

 figure(2)
 hold on
%plot(r,p_2(1:length(r)), 'r-')                          % Analítico Akhnoukh com fonte pontual
plot(r,p_2_J(1:length(r)), 'k-')                        % Analítico Jacobsen com fonte pontual
  plot(x,p_flow, 'g-')                                  % Analítico Akhnoukh com escoamento
% plot(Simulation(:,8), (Simulation(:,5))-101325, 'r-')   % Simulação
hold off
 xlabel('x (m)')
ylabel('p (Pa)')