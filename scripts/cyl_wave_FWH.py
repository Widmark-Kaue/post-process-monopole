#%% Librarys
from time import time
from unittest import skip
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel2
from pathlib import Path
#%% FUNÇÕES
def importData(simulation:str, probe:int=2, time:float = None)-> tuple:
    PATH_DATA   = Path().absolute().parent / 'data'
    PROBES      = PATH_DATA / 'probes' / simulation /'p.txt'
    FWH         = PATH_DATA / 'pressureData' / simulation / 'FWH_surfacePressure.txt'

    if time == None:
        t, p         = np.loadtxt(PROBES, usecols = (0, probe+1), unpack=True)
        fwh_t, fwh_p = np.loadtxt(FWH, usecols=(0,probe+1), unpack=True)
        print(f'Probe: {probe}')
    else:
        t   = np.loadtxt(PROBES / simulation /'p.txt', usecols = 0)
        pos = t.searchsorted(time)
        p   = np.loadtxt(PROBES, skiprows=len(t)-pos) 
        print(f'Time = {t[pos]}')
    return ((t,p), (fwh_t, fwh_p)) 

def plot(FWH:tuple, SIM:tuple, fanalitic:any) -> None:
    fwh_t, fwh_p    = FWH
    t, p            = SIM
   
    if fanalitic != None:
        p_2_E           = fanalitic(t) 
        plt.plot(t, p_2_E, 'k--', label = 'Analítico')
        plt.plot(t, p - 101325,'r--', label = 'Cálculo Direto', alpha = 0.5)
    else:
        plt.plot(t, p - 101325,'r--', label = 'FWH2')
    plt.plot(fwh_t, fwh_p - 101325, 'b--', label = 'FWH', alpha = 0.35)

    plt.xlabel('Tempo [s]')
    plt.ylabel('Pressão [Pa]')
    plt.legend()
    
    plt.grid()
    plt.show()
#%% Constantes
microphone = np.arange(2, 112, 10)
#%% Solução Analítica
rf = 0.05715 / 2
S = 0.1
c0 = 340.29
c = 331.45
freq = 100
omega = freq * 2 * np.pi
T0 = 273.15
T = (c0 / c) ** 2 * T0
rho0 = 101325 / (287.058 * T) #% Cálculo da densidade
area = 2 * np.pi * rf
velocity = S / area

H_1_fonte_J = hankel2(1, (omega * rf / c0))
A0 = velocity*1j*rho0*c0/H_1_fonte_J
H_0_E = lambda x: hankel2(0, omega*x/c0)
P_2_E = lambda x: A0*H_0_E(x)
p_2_E = lambda t, x=microphone[2]: P_2_E(x)*np.exp(1j*omega*t).imag

#%% Importando dados
# # SIM, FWH40 = importData('backward40')
# _,FWH20    = importData('backward20')
SIM,FWH2     = importData('monopole2')

#%% Plot
plot(FWH2, SIM, fanalitic = p_2_E)

# %% testes
# rexp = np.linspace(rvec[0],rvec[-1])
# plt.plot(rexp, p_2_E(t = 0.5,x = rexp), 'k', label = 'analítico')

# plt.xlabel('r [m]')
# plt.ylabel('Pressão [Pa]')
# plt.legend()
# plt.grid()
# plt.show()
# %%
