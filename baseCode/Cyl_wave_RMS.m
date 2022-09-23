% clear all
close all
clc

% save('pLimited16ppw.txt','pLimited16ppw','-ascii')
% save('p32ppw.txt','p32ppw','-ascii')

%% ERROS RMS

r = 2:10:102;
%% ERRO REFINO

  Simulation8 = load('p8ppw.txt');       % resultado da simulação
%  Simulation16 = load('p16ppw.txt');
% Simulation32 = load('p32ppw.txt');
% Simulation64 = load('p64ppw.txt');
% 
% 
  RMS_8ppw = rms(Simulation8,r)*100;
%  RMS_16ppw = rms(Simulation16,r)*100;
% RMS_32ppw = rms(Simulation32,r)*100;
% RMS_64ppw = rms(Simulation64,r)*100;
% 
% RMS = [RMS_8ppw, RMS_16ppw, RMS_32ppw, RMS_64ppw];
% PPW = [8, 16, 32, 64];
% 
% figure(2)
% hold on
% plot(PPW,RMS,'-o','color',[0.6350 0.0780 0.1840])
% hold off
% xticks([8 16 32 64]);
% xlabel('PPW')
% ylabel('RMS (%)')


%% ERRO TEMPORAL

% pBackward10 = load('pBackward10.txt');
% pBackward20 = load('pBackward20.txt');
% pBackward40 = load('pBackward40.txt');
% pBackward80 = load('pBackward80.txt');
% pBackward4000 = load('p32ppw.txt');
% RMS_Backward10 = rms(pBackward10,r)*100;
% RMS_Backward20 = rms(pBackward20,r)*100;
% RMS_Backward40 = rms(pBackward40,r)*100;
% RMS_Backward80 = rms(pBackward80,r)*100;
% RMS_Backward4000 = rms(pBackward4000,r)*100;
% 
% RMS_Backward = [RMS_Backward10, RMS_Backward20, RMS_Backward40, RMS_Backward80];
% 
% pEuler10 = load('pEuler10.txt');
% pEuler20 = load('pEuler20.txt');
% pEuler40 = load('pEuler40.txt');
% pEuler80 = load('pEuler80.txt');
% pEuler4000 = load('pEuler.txt');
% RMS_Euler10 = rms(pEuler10,r)*100;
% RMS_Euler20 = rms(pEuler20,r)*100;
% RMS_Euler40 = rms(pEuler40,r)*100;
% RMS_Euler80 = rms(pEuler80,r)*100;
% RMS_Euler4000 = rms(pEuler4000,r)*100;
% 
% RMS_Euler = [RMS_Euler10, RMS_Euler20, RMS_Euler40, RMS_Euler80];
% 
% n = [10, 20, 40, 80];
% 
% figure
% hold on
% plot(n,RMS_Euler,'-o','color','k')
% plot(n,RMS_Backward,'--o','color','b')
% hold off
% xticks([10 20 40 80]);
% xlabel('n')
% ylabel('RMS (%)')
% legend('Euler', 'Backward', 'Location','best');


%% ERRO ESQUEMAS NUMERICOS

% pLimited16ppw = load('pLimited16ppw.txt');        % resultado da simulação
% pLimited32ppw = load('pLimited32ppw.txt');
% pLinearUpwind16ppw = load('pLinearUpwind16ppw.txt');
% pLinearUpwind32ppw = load('pLinearUpwind32ppw.txt');
% pUpwind16ppw = load('pUpwind16ppw.txt');
% pUpwind32ppw = load('pUpwind32ppw.txt');
% pVanLeer16ppw = load('p16ppw.txt');
% pVanLeer32ppw = load('p32ppw.txt');
% 
% RMS_limitedLinear16 = rms(pLimited16ppw,r)*100;
% RMS_limitedLinear32 = rms(pLimited32ppw,r)*100;
% RMS_linearUpwind16 = rms(pLinearUpwind16ppw,r)*100;
% RMS_linearUpwind32 = rms(pLinearUpwind32ppw,r)*100;
% RMS_Upwind16 = rms(pUpwind16ppw,r)*100;
% RMS_Upwind32 = rms(pUpwind32ppw,r)*100;
% RMS_vanLeer16 = rms(pVanLeer16ppw,r)*100;
% RMS_vanLeer32 = rms(pVanLeer32ppw,r)*100;
% 
% 
% %% ESQUEMA X RMS
% 
% RMS_esquema16 = [RMS_limitedLinear16, RMS_linearUpwind16, RMS_Upwind16, RMS_vanLeer16];
% RMS_esquema32 = [RMS_limitedLinear32, RMS_linearUpwind32, RMS_Upwind32, RMS_vanLeer32];
% 
% % RMS_limitedLinear = [RMS_limitedLinear16, RMS_limitedLinear32];
% % RMS_linearUpwind = [RMS_linearUpwind16, RMS_linearUpwind32];
% % RMS_Upwind = [RMS_Upwind16, RMS_Upwind32];
% % RMS_vanLeer = [RMS_vanLeer16, RMS_vanLeer32];
% % 
% % esquemas_PPW = [16, 32];
% 
% esquemas_label=["limitedLinear","linearUpwind","upwind","vanLeer"];
% esquemas = categorical(esquemas_label);
% 
% figure(5)
% hold on
% plot(esquemas,RMS_esquema16,'o','color','r')
% plot(esquemas,RMS_esquema32,'*','color','k')
% hold off
% grid on
% xlabel('Esquemas Numericos')
% ylabel('RMS (%)')
% legend('16 PPW', '32 PPW', 'Location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(7)
% hold on
% plot(pBackward20(:,1),(pBackward20(:,11)-101325),'-','color','r')
% hold off
% xlabel('t (s)')
% ylabel('p (Pa)')

function RMS_max = rms(Simulation,r)
    num_probes = size(Simulation);
    num_probes = num_probes(1,2)-1;
    RMS_total = zeros(num_probes,1);
    for i=1:num_probes
        RMS_total(i) = rms_error(Simulation,r(i),i+1);
    end
    RMS_max = max(RMS_total);
end

function RMS = rms_error(Simulation,r,probe_column)
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
    t_transiente = r/c_0 + 5/freq;
    
    count = 0;
    for i=1:length(Simulation(:,1))
        if Simulation(i,1) >= t_transiente
            count = count+1;
            t_erro(count)=Simulation(i,1);
            p_sim(count) = Simulation(i,probe_column);
        end
    end

    
    H_0_E = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
    P_2_E = A0*H_0_E;
    p_2_E = imag(P_2_E*exp(1i*omega*t_erro));
    
    figure
    hold on
    plot(t_erro,p_2_E,'-','color','r')
    plot(t_erro,(p_sim-101325),'-','color','k')
    hold off
    xlabel('t (s)')
    ylabel('p (Pa)')
    legend('Analítico', 'Direto', 'Location','best');
                       
    % RMS
    Num1 = (p_sim-101325)-p_2_E;
    Num2 = Num1.^2;
    Num = trapz(t_erro,Num2);
    Den = trapz(t_erro,p_2_E.^2);
    RMS = Num/Den;
end 