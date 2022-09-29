#%% Librarys
from cProfile import label
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import hankel2
from scipy.integrate import trapz

#%% FUNÇÕES
def pressure(r:float or np.ndarray = None, 
            t:float or np.ndarray = None) -> float or np.ndarray:
    rf          = 0.05715 / 2
    S           = 0.1
    c0          = 340.29
    c           = 331.45
    freq        = 100
    omega       = freq * 2 * np.pi
    T0          = 273.15
    T           = (c0 / c) ** 2 * T0
    rho0        = 101325 / (287.058 * T) #% Eq. dos Gases ideias
    area        = 2 * np.pi * rf
    velocity    = S / area

    H_1_fonte_J = hankel2(1, (omega * rf / c0))
    A0          = velocity*1j*rho0*c0/H_1_fonte_J
    
    H_0_E       = lambda r1: hankel2(0, omega*r1/c0)
    P_2_E       = lambda r1: A0*H_0_E(r1)
    p_2_E       = lambda t1,r1: P_2_E(r1)*np.exp(1j*omega*t1).imag
    
    if r == None and t == None:
        return p_2_E
    elif any(i != None for i in [t,r]):
        aux     = t == None
        aux2    = lambda t1: p_2_E(t1,r1=r) if aux else lambda r1: p_2_E(r1, t1=t)
        return aux2
    else:
        return p_2_E(t,r)



def importData(simulation:str, probe:int=2, time:float = None)-> tuple:
    PATH_DATA   = Path().absolute().parent / 'data'
    PROBES      = PATH_DATA / 'probes' / simulation /'p.txt'
    FWH         = PATH_DATA / 'pressureData' / simulation / 'FWH_surfacePressure.txt'

    toPa        = 101325

    if time == None:
        t, p         = np.loadtxt(PROBES, usecols = (0, probe+1), unpack=True)
        fwh_t, fwh_p = np.loadtxt(FWH, usecols=(0,probe+1), unpack=True)
        print(f'Probe: {probe}')
        return ((t,p-toPa), (fwh_t, fwh_p-toPa)) 
    else:
        t   = np.loadtxt(PROBES / simulation /'p.txt', usecols = 0)
        pos = t.searchsorted(time)
        p   = np.loadtxt(PROBES, skiprows=pos-1)[1:] 
        print(f'Time = {t[pos]}')
        return p-toPa


def rms(tp:tuple, pfunc:any, r:int, 
        freq:float = 100, c0:float = 340.29, plot: bool = True):
        
    
    t, p                = tp
    transientTime       = r/c0 +5/freq
    fil                 = lambda x: x >= transientTime
    efectiveTime        = np.array(list(filter(fil, t)))
    efectivePressure    = np.array([p[t.index(i)] for i in efectiveTime]) 

    #% Plot
    if plot:
        plt.plot(efectiveTime, pfunc(efectiveTime), 'k-', label = 'Analítico')
        plt.plot(efectiveTime, efectivePressure, 'r-', label = 'Direto')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Pressão [Pa]')
        plt.legend()
        plt.grid()
        plt.show()
    
    #% RMS
    aux1    = (p - pfunc(efectiveTime))**2
    rms     = trapz(aux1, efectiveTime)/trapz(pfunc(efectiveTime), efectiveTime)
    
    return rms  


def plotTime(FWH:tuple, SIM:tuple, fanalitic:any) -> None:
    fwh_t, fwh_p    = FWH
    t, p            = SIM
   
    if fanalitic != None:
        p_2_E           = fanalitic(t) 
        plt.plot(t, p_2_E, 'k', label = 'Analítico')
        plt.plot(t, p,'r-.', label = 'Cálculo Direto', alpha = 0.5)
    else:
        plt.plot(t, p,'r-.', label = 'FWH2')
    plt.plot(fwh_t, fwh_p, 'b--', label = 'FWH', alpha = 0.35)

    plt.xlabel('Tempo [s]')
    plt.ylabel('Pressão [Pa]')
    plt.legend()
    
    plt.grid()
    plt.show()
#%% Parâmetros de simulação
microphone = np.arange(2, 112, 10)

#%% Importando dados
# # SIM, FWH40 = importData('backward40')
# _,FWH20    = importData('backward20')
SIM,FWH2     = importData('monopole2')

#%% Plot
p_time = lambda t: pressure(r = microphone[2], t=t) 
plotTime(FWH2, SIM, fanalitic = p_time)

# %% testes
# rexp = np.linspace(rvec[0],rvec[-1])
# plt.plot(rexp, p_2_E(t = 0.5,x = rexp), 'k', label = 'analítico')

# plt.xlabel('r [m]')
# plt.ylabel('Pressão [Pa]')
# plt.legend()
# plt.grid()
# plt.show()
# %%
