% clear all
close all
clc


% save('pUpwind32ppw.txt','pUpwind32ppw','-ascii','-double')

%% ERROS RMS

% tipo = 'temp';
% tipo = 'refino';
 tipo = 'esquemas';
% tipo = 0;

r = 2:10:102;
%% ERRO REFINO

if strcmp(tipo,'refino') == 1

% Simulation8 = load('p8ppw.txt');      
Simulation16 = load('p16ppw.txt');
Simulation32 = load('p32ppw.txt');
Simulation64 = load('p64ppw.txt');


%  RMS_8ppw = dissip8(Simulation8,r)*100;
RMS_16ppw = dissip(Simulation16,r)*100;
RMS_32ppw = dissip(Simulation32,r)*100;
RMS_64ppw = dissip(Simulation64,r)*100;

RMS = [RMS_16ppw, RMS_32ppw, RMS_64ppw];

p_refino = plot_dissip(RMS, RMS, tipo);

end

%% ERRO TEMPORAL

if strcmp(tipo,'temp') == 1

pBackward10 = load('pBackward10.txt');
pBackward20 = load('pBackward20.txt');
pBackward40 = load('pBackward40.txt');
pBackward80 = load('pBackward80.txt');
pBackward4000 = load('p32ppw.txt');
RMS_Backward10 = dissip(pBackward10,r)*100;
RMS_Backward20 = dissip(pBackward20,r)*100;
RMS_Backward40 = dissip(pBackward40,r)*100;
RMS_Backward80 = dissip(pBackward80,r)*100;
RMS_Backward4000 = dissip(pBackward4000,r)*100;

RMS_Backward = [RMS_Backward10, RMS_Backward20, RMS_Backward40, RMS_Backward80];

 pEuler10 = load('pEuler10.txt');
pEuler20 = load('pEuler20.txt');
pEuler40 = load('pEuler40.txt');
pEuler80 = load('pEuler80.txt');
pEuler4000 = load('pEuler.txt');
 RMS_Euler10 = dissip(pEuler10,r)*100;
RMS_Euler20 = dissip(pEuler20,r)*100;
RMS_Euler40 = dissip(pEuler40,r)*100;
RMS_Euler80 = dissip(pEuler80,r)*100;
RMS_Euler4000 = dissip(pEuler4000,r)*100;

RMS_Euler = [RMS_Euler10, RMS_Euler20, RMS_Euler40, RMS_Euler80];

p_temp = plot_dissip(RMS_Euler, RMS_Backward, tipo);

end


%% ERRO ESQUEMAS NUMERICOS

if strcmp(tipo,'esquemas') == 1

pLimited16ppw = load('pLimited16ppw.txt');        % resultado da simulação
pLimited32ppw = load('pLimited32ppw.txt');
pLinearUpwind16ppw = load('pLinearUpwind16ppw.txt');
pLinearUpwind32ppw = load('pLinearUpwind32ppw.txt');
pUpwind16ppw = load('pUpwind16ppw.txt');
pUpwind32ppw = load('pUpwind32ppw.txt');
pVanLeer16ppw = load('p16ppw.txt');
pVanLeer32ppw = load('p32ppw.txt');

RMS_limitedLinear16 = dissip(pLimited16ppw,r)*100;
RMS_limitedLinear32 = dissip(pLimited32ppw,r)*100;
RMS_linearUpwind16 = dissip(pLinearUpwind16ppw,r)*100;
RMS_linearUpwind32 = dissip(pLinearUpwind32ppw,r)*100;
RMS_Upwind16 = dissip(pUpwind16ppw,r)*100;
RMS_Upwind32 = dissip(pUpwind32ppw,r)*100;
RMS_vanLeer16 = dissip(pVanLeer16ppw,r)*100;
RMS_vanLeer32 = dissip(pVanLeer32ppw,r)*100;


%% ESQUEMA X RMS

RMS_esquema16 = [RMS_limitedLinear16, RMS_linearUpwind16, RMS_Upwind16, RMS_vanLeer16];
RMS_esquema32 = [RMS_limitedLinear32, RMS_linearUpwind32, RMS_Upwind32, RMS_vanLeer32];

p_temp = plot_dissip(RMS_esquema16, RMS_esquema32, tipo);

end



function RMS_total = plot_dissip(Vector1, Vector2, tipo)
    num_probes1 = size(Vector1);
    num_probes = num_probes1(1,1);
    RMS_total = zeros(num_probes,1);
    for i=1:num_probes
        % REFINO
        if strcmp(tipo,'refino') == 1
            PPW = [16, 32, 64];

            figure
            hold on
            plot(PPW,Vector1(i,:),'-o','color',[0.6350 0.0780 0.1840])
            hold off
            xticks([16 32 64]);
            xlabel('PPW')
            ylabel('Erro de amplitude (%)')
            set(gca,'FontSize',20)
        end
        % TEMPORAL
        if strcmp(tipo,'temp') == 1
            
            n = [10, 20, 40, 80];

            figure
            hold on
            plot(n,Vector1(i,:),'-o','color','k')
            plot(n,Vector2(i,:),'--o','color','b')
            hold off
            xticks([10 20 40 80]);
            xlabel('n')
            ylabel('Erro de amplitude (%)')
            legend('Euler', 'Backward', 'Location','best');
            set(gca,'FontSize',20)
        end
        % ESQUEMAS
        if strcmp(tipo,'esquemas') == 1
            esquemas_label=["limitedLinear","linearUpwind","upwind","vanLeer"];
            esquemas = categorical(esquemas_label);

            figure
            hold on
            bar(esquemas,Vector1(i,:),'c')
            bar(esquemas,Vector2(i,:),'r')
            hold off
            grid on
            xlabel('Esquemas Numericos')
            ylabel('Erro de amplitude (%)')
            legend('16 PPW', '32 PPW', 'Location','best');
            set(gca,'FontSize',20)
        end
    end
end

function RMS_total = dissip(Simulation,r)
    num_probes1 = size(Simulation);
    num_probes = num_probes1(1,2)-1;
    RMS_total = zeros(num_probes,1);
    for i=1:num_probes
        RMS_total(i) = dissip_probe(Simulation,r(i),i+1);
    end
    %RMS_max = max(RMS_total);
end

function Erro_diss = dissip_probe(Simulation,r,probe_column)
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
    
     tempos = length(Simulation(:,1));
%     pico = max(Simulation(:,probe_column));
%     
%     % Analítico
%     for i=1:tempos
%         if Simulation(i,probe_column) == pico
%             t_transiente = Simulation(i,1);
%             break
%         end
%     end
    
    t_transiente = r/c_0 + 5/freq;
    count = 0;
    for i=1:tempos
        if Simulation(i,1) >= t_transiente
            count = count+1;
            t_erro(count) = Simulation(i,1);
            p_sim(count) = Simulation(i,probe_column);
        end
    end

    
    H_0_E = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
    P_2_E = A0*H_0_E;
    p_2_E = imag(P_2_E*exp(1i*omega*t_erro));
    
%     figure
%     hold on
%     plot(t_erro,p_2_E,'-','color','r')
%     plot(t_erro,(p_sim-101325),'-','color','k')
%     hold off
%     xlabel('t (s)')
%     ylabel('p (Pa)')
%     legend('Analítico', 'Direto', 'Location','best');
                       
    % Erro
    prms_sim = trapz(t_erro,(p_sim-101325).^2);
    prms_anal = trapz(t_erro,p_2_E.^2);
    Erro_diss = abs(prms_sim-prms_anal)/prms_anal;
end 

function RMS_total = dissip8(Simulation,r)
    num_probes1 = size(Simulation);
    num_probes = num_probes1(1,2)-1;
    RMS_total = zeros(num_probes,1);
    for i=1:num_probes
        RMS_total(i) = dissip_probe8(Simulation,r(i),i+1);
    end
    %RMS_max = max(RMS_total);
end

function Erro_diss = dissip_probe8(Simulation,r,probe_column)
    if probe_column < 12
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
    
     tempos = length(Simulation(:,1));
%     pico = max(Simulation(:,probe_column));
%     
%     % Analítico
%     for i=1:tempos
%         if Simulation(i,probe_column) == pico
%             t_transiente = Simulation(i,1);
%             break
%         end
%     end
    
    [pico_sim1, pos_max_sim] = max(Simulation(:,probe_column));
     t_transiente = (Simulation(pos_max_sim,1));

    %t_transiente = r/c_0 + 5/freq;
    count = 0;
    for i=1:tempos
        if Simulation(i,1) >= t_transiente
            count = count+1;
            t_erro(count) = Simulation(i,1);
            p_sim(count) = Simulation(i,probe_column);
        end
    end

    
    H_0_E = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
    P_2_E = A0*H_0_E;
    p_2_E = imag(P_2_E*exp(1i*omega*t_erro));
    
%     figure
%     hold on
%     plot(t_erro,p_2_E,'-','color','r')
%     plot(t_erro,(p_sim-101325),'-','color','k')
%     hold off
%     xlabel('t (s)')
%     ylabel('p (Pa)')
%     legend('Analítico', 'Direto', 'Location','best');
                       
    % Erro
    prms_sim = trapz(t_erro,(p_sim-101325).^2);
    prms_anal = trapz(t_erro,p_2_E.^2);
    Erro_diss = abs(prms_sim-prms_anal)/prms_anal;
    else 
         Erro_diss = NaN;
    end
end 