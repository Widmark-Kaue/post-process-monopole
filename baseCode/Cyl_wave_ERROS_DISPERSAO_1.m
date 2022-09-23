% clear all
close all
clc

% save('pLimited16ppw.txt','pLimited16ppw','-ascii')

%% ERROS RMS

% tipo = 'temp';
% tipo = 'refino';
 tipo = 'esquemas';

r = 2:10:102;
%% ERRO REFINO

if strcmp(tipo,'refino') == 1

%Simulation8 = load('p8ppw.txt');      
Simulation16 = load('p16ppw.txt');
 Simulation32 = load('p32ppw.txt');
%%%Simulation32 = load('pBackward80.txt');
Simulation64 = load('p64ppw.txt');

%  disp_probe(Simulation64,102,12);
% disp(Simulation32,r);


% RMS_8ppw = disp8(Simulation8,r);
RMS_16ppw = disp(Simulation16,r);
RMS_32ppw = disp(Simulation32,r);
RMS_64ppw = disp(Simulation64,r);

RMS = [RMS_16ppw, RMS_32ppw, RMS_64ppw];

p_refino = plot_disp(RMS, RMS, tipo);

end

%% ERRO TEMPORAL

if strcmp(tipo,'temp') == 1

pBackward10 = load('pBackward10.txt');
 pBackward20 = load('pBackward20.txt');
pBackward40 = load('pBackward40.txt');
pBackward80 = load('pBackward80.txt');
pBackward4000 = load('p32ppw.txt');
RMS_Backward10 = disp(pBackward10,r);
 RMS_Backward20 = disp(pBackward20,r);
RMS_Backward40 = disp(pBackward40,r);
RMS_Backward80 = disp(pBackward80,r);
RMS_Backward4000 = disp(pBackward4000,r);

RMS_Backward = [RMS_Backward10, RMS_Backward20, RMS_Backward40, RMS_Backward80];

pEuler10 = load('pEuler10.txt');
pEuler20 = load('pEuler20.txt');
pEuler40 = load('pEuler40.txt');
  pEuler80 = load('pEuler80.txt');
pEuler4000 = load('pEuler.txt');
RMS_Euler10 = disp(pEuler10,r);
RMS_Euler20 = disp(pEuler20,r);
RMS_Euler40 = disp(pEuler40,r);
  RMS_Euler80 = disp(pEuler80,r);
RMS_Euler4000 = disp(pEuler4000,r);

RMS_Euler = [RMS_Euler10, RMS_Euler20, RMS_Euler40, RMS_Euler80];

p_temp = plot_disp(RMS_Euler, RMS_Backward, tipo);

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

RMS_limitedLinear16 = disp(pLimited16ppw,r);
RMS_limitedLinear32 = disp(pLimited32ppw,r);
RMS_linearUpwind16 = disp(pLinearUpwind16ppw,r);
RMS_linearUpwind32 = disp(pLinearUpwind32ppw,r);
RMS_Upwind16 = disp(pUpwind16ppw,r);
RMS_Upwind32 = disp(pUpwind32ppw,r);
RMS_vanLeer16 = disp(pVanLeer16ppw,r);
RMS_vanLeer32 = disp(pVanLeer32ppw,r);


%% ESQUEMA X RMS

RMS_esquema16 = [RMS_limitedLinear16, RMS_linearUpwind16, RMS_Upwind16, RMS_vanLeer16];
RMS_esquema32 = [RMS_limitedLinear32, RMS_linearUpwind32, RMS_Upwind32, RMS_vanLeer32];

p_temp = plot_disp(RMS_esquema16, RMS_esquema32, tipo);

end



function RMS_total = plot_disp(Vector1, Vector2, tipo)
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
            ylabel('Erro de fase (grau)')
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
            ylabel('Erro de fase (grau)')
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
            ylabel('Erro de fase (grau)')
            legend('16 PPW', '32 PPW', 'Location','best');
            set(gca,'FontSize',20)
        end
    end
end

function RMS_total = disp(Simulation,r)
    num_probes1 = size(Simulation);
    num_probes = num_probes1(1,2)-1;
    RMS_total = zeros(num_probes,1);
    for i=1:num_probes
        RMS_total(i) = disp_probe(Simulation,r(i),i+1);
    end
    %RMS_max = max(RMS_total);
end

function Erro_dispersao = disp_probe(Simulation,r,probe_column)
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
    
     tempos_sim = length(Simulation(:,1));
    passo = find(Simulation(:,1)==0.01);
%     
%     % Analítico
%     
    
    t_transiente = r/c_0 + 5/freq;
    count = 0;
    for i=1:tempos_sim
        if Simulation(i,1) >= t_transiente
            count = count+1;
            t_erro(count) = Simulation(i,1);
            p_sim(count) = Simulation(i,probe_column);
        end
    end

    
    H_0_E = besselh(0, 2, (omega*r/c_0));                  %Função de Hankel
    P_2_E = A0*H_0_E;
    p_2_E = imag(P_2_E*exp(1i*omega*t_erro));
    
%     for i=1:length(t_erro)
%         if p_2_E(i) == max(p_2_E)
%             t_analitico = t_erro(i);
%             break
%         end
%     end
%      pos_analitico = find(abs(Simulation(:,1)-t_analitico) < 0.000001); 
% 

    for i=passo:length(t_erro)
        if p_2_E(i)< 0.000001
            t0 = t_erro(i);
            pos_0 = i;
%             for j = i:-1:1
%                 if p_2_E(j) == max(p_2_E)
%                 t_analitico = t_erro(j);
%                 break
%                 end
%             end
            break
        end
    end
%     
%     
%      pos_analitico = find(abs(Simulation(:,1)-t_analitico) < 0.000001); 
 %    pos_0 = find(abs(Simulation(:,1)-t0) < 0.000001);
     
     
     for i=(pos_0):-1:(pos_0-passo)
         vetor_anal(i) = p_2_E(i);
     end
    
    
     pico_anal = max(vetor_anal);
     
    
    
    for i=(pos_0):-1:(pos_0-passo)
            if p_2_E(i) == pico_anal
                t_analitico = t_erro(i);
                break
            end
    end

     
     for i=(pos_0):-1:(pos_0-passo)
         vetor_pico(i) = p_sim(i);
     end
    
    
     pico_sim = max(vetor_pico);
     
    
    
    for i=(pos_0):-1:(pos_0-passo)
            if p_sim(i) == pico_sim
                t_pico_sim2 = t_erro(i);
                break
            end
    end
    
    
%     figure
%     hold on
%     plot(t_erro,p_2_E,'-','color','r')
%     plot(t_erro,(p_sim-101325),'-','color','k')
%     hold off
%     xlabel('t (s)')
%     ylabel('p (Pa)')
%     legend('Analítico', 'Direto', 'Location','best');
%     set(gca,'FontSize',20)
    
    Deltat2 = abs(t_analitico - t_pico_sim2);
    Erro_dispersao =  rad2deg(Deltat2*omega);
end 

