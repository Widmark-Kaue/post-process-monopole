#%% Library
import numpy as np
import wave as wv
import matplotlib.pyplot as plt


#%% Teste 3

# Comportamento espacial
pSim,pFWH,pFWH2 = wv.importData('pimpleT3', time = 0.5, case = 'monopoleFlow')
x = np.linspace(2,102,len(pSim))

plt.plot(x,pSim, label = 'Cálculo Direto')
plt.xlabel('x [m]')
plt.ylabel('p [Pa]')
plt.title('Teste 3\nTempo 0.5 s')
plt.grid()
plt.legend()
plt.show()

#%% Comportamento temporal
microphone = np.linspace(2,102, 11)
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('pimpleT3', case = 'monopoleFlow', probe=m)

    ti = f'Teste 3\nComportamento temporal (r = {round(microphone[m],2)} [m])'
   
    plt.plot(SIM[0], SIM[1], 'r--', label = 'Cálculo direto', alpha = 0.5)
    plt.plot(FWH[0], FWH[1], 'b--', label = 'FWH1', alpha = 0.85)
    plt.plot(FWH2[0], FWH[1], 'g--', label = 'FWH2', alpha = 0.75)

    plt.xlabel('t [s]')
    plt.ylabel('P [Pa]') 
    plt.title(ti)

    plt.legend()
    plt.grid()
    plt.show()

#%% Teste 6
microphone = np.linspace(2,102, 11)
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('pimpleT6', case = 'monopoleFlow', probe=m)

    ti = f'Teste6\nComportamento temporal (r = {round(microphone[m],2)} [m])'
   
    plt.plot(SIM[0], SIM[1], 'r--', label = 'Cálculo direto', alpha = 0.5)
    plt.plot(FWH[0], FWH[1], 'b--', label = 'FWH1', alpha = 0.85)
    plt.plot(FWH2[0], FWH[1], 'g--', label = 'FWH2', alpha = 0.75)

    plt.xlabel('t [s]')
    plt.ylabel('P [Pa]') 
    plt.title(ti)

    plt.legend()
    plt.grid()
    plt.show()


#%% teste função
x = np.linspace(-102,102,10000)
pFlow = wv.pressureFlow(x = x, t = 0.5, M = 0.1)

plt.plot(x,pFlow, label = "analytical")
plt.xlabel("x [m]")
plt.ylabel("p [Pa]")
plt.title("t = 0.5")
plt.grid()
plt.legend()
plt.plot()

# %% Teste 6
