% clear all
close all
clc

%save('pBackward40.txt','pBackward40','-ascii','-double')
%save('FWH1Backward40.txt','FWHtime1','-ascii','-double')
%save('FWH2Backward40d-1.txt','FWH2time1','-ascii','-double')
%type('FWH1Backward20.txt')

%% TCC - FWH

r = 2:10:102; % coordenadas das probes

% pBackward20 = load('pBackward20.txt');
% FWH1Backward20 = load('FWH1Backward20.txt');
% FWH2Backward20 = load('FWH2Backward20.txt');

pBackward40 = load('pBackward40.txt');
FWH1Backward40 = load('FWH1Backward40.txt');
FWH2Backward40 = load('FWH2Backward40d-1.txt');


% Backward20FWH = plot_comp(pBackward20,FWH2Backward20,r,7);
Backward40FWH = plot_comp(pBackward40,FWH2Backward40,r,2);
%ppw16FWH = plot_comp(ps01ppw16,FWHtimes01ppw16,r,1);


function FWH = plot_comp(Simulation,FWH,r,probe)
    probe_column = probe+1;
    r = r(probe);
    r_fonte = 0.05715/2; % raio da fonte
    S = 0.1;          % amplitude da vazão definida na simulação
    c_0=340.29;       % velocidade do som utilizada no cálculo do comprimento de onda
    freq=100;         % frequência da onda 
    omega = freq*2*pi;
    c = 331.45;       % velocidade do som a 0 ºC
    T0 = 273.15;      % 0 ºC em Kelvin
    T = (c_0/c)^2*T0; % temperatura da simulação 

    % CÁLCULO DA DENSIDADE
    rho0 = 101325/(287.058*T); % equação dos gases ideais
    
    Area = 2*pi*r_fonte;
    Velocity = S/Area;
    H_1_fonte_J = besselh(1, 2, (omega*r_fonte/c_0));

    A0 = Velocity*1i*rho0*c_0/H_1_fonte_J;
    
    % Analítico
    t = Simulation(:,1);
    H_0_E = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
    P_2_E = A0*H_0_E;
    p_2_E = imag(P_2_E*exp(1i*omega*t));
           
    figure
    hold on
    plot(t,p_2_E(1:length(t)), 'k-')                        % Analítico Jacobsen com fonte pontual
    plot(t, (Simulation(:,probe_column))-101325, 'r-')   % Simulação
    plot(FWH(:,1), (FWH(:,probe_column)), 'b--')   % Simulação
    hold off
    xlabel('t (s)')
    ylabel('p (Pa)')
    legend('Analítico', 'Cálculo Direto', 'Analogia Acústica', 'Location','best');
end 