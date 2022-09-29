#%% Librarys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.special import hankel2

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
    p_2_E       = lambda t1,r1: (P_2_E(r1)*np.exp(1j*omega*t1)).imag
    
    if r == None and t == None:
        return p_2_E
    elif any(i != None for i in [t,r]):
        aux     = t == None
        aux2    = lambda t1: p_2_E(t1,r1=r) if aux else lambda r1: p_2_E(r=r1, t1=t)
        return aux2
    else:
        return p_2_E(t,r)



def importData(simulation:str, probe:int=2, time:float = None)-> tuple:
    
    assert 0<= probe <= 10, 'ERRO: O valor da probe deve estar entre 0 e 10'
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
        t   = np.loadtxt(PROBES, usecols = 0)
        pos = np.searchsorted(t,time)-1
        p   = np.loadtxt(PROBES, skiprows=12 + pos)[1,1:] 
        print(f'Time = {t[pos]} \nPos = {pos}')
        return p - toPa


def rmsTime(tp:tuple, robs:float, freq:float = 100, 
            c0:float = 340.29, plot: bool = True)-> float:
        
    
    t, p                = tp
    transientTime       = robs/c0 +5/freq
    fil                 = lambda x: x >= transientTime
    efectiveTime        = np.array(list(filter(fil, t)))
    t0                  = list(t).index(efectiveTime[0])
    efectivePressure    = p[t0:] 

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
    pfunc   = pressure(r = robs)
    aux1    = (efectivePressure - pfunc(efectiveTime))**2
    num     = trapz(aux1, efectiveTime)
    den     = trapz(pfunc(efectiveTime)**2, efectiveTime)
    rms     = num/den

    return rms  

def rmsSpacial(rp:tuple, tobs:float = 0.5, plot:bool = True):
    
    r, p    = rp
    pfunc   = pressure(t = tobs)
    
    if plot:
        plt.plot(r, pfunc(r), 'k-', label = 'Analítico')
        plt.plot(r, p, 'r-', label = 'Direto')
        plt.xlabel('raio [m]')
        plt.ylabel('Pressão [Pa]')
        plt.legend()
        plt.grid()
        plt.show()
    
    aux1    = (p - pfunc(r))**2
    num     = trapz(aux1, r)
    den     = trapz(pfunc(r)**2, r)
    rms     = num/den 

    return rms

def plotTime(FWH:tuple, SIM:tuple, fanalitic:any, title:str = None) -> None:
    fwh_t, fwh_p    = FWH
    t, p            = SIM
   
    if fanalitic != None:
        p_2_E           = fanalitic(t) 
        plt.plot(t, p_2_E, 'k', label = 'Analítico')
        plt.plot(t, p,'r-.', label = 'Cálculo Direto', alpha = 0.5)
    else:
        plt.plot(t, p,'r-.', label = 'FWH2')
    plt.plot(fwh_t, fwh_p, 'b--', label = 'FWH', alpha = 0.35)

    if title!=None: plt.title(title)
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
SIM,FWH2     = importData('monopole2', probe = 10)
pr           = importData('monopole2', time  = 0.5)

#%% Plot
p_time = pressure(r  = microphone[2])
plotTime(FWH2, SIM, fanalitic = p_time)


#%% RMS
RMS = rmsTime(SIM, r = microphone[2])*100
print(RMS)

# %% teste
a = np.loadtxt('teste.txt',skiprows=3)
print(a)
# %%
