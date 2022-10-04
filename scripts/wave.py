#%% Librarys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.special import hankel2

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
    Lambda      = c/freq

    H_1_fonte_J = hankel2(1, (omega * rf / c0))
    A0          = velocity*1j*rho0*c0/H_1_fonte_J
    
    H_0_E       = lambda r1: hankel2(0, omega*r1/c0)
    P_2_E       = lambda r1: A0*H_0_E(r1)
    p_2_E       = lambda t1,r1: (P_2_E(r1)*np.exp(1j*omega*t1)).imag
    
    if r == None and t == None:
        return p_2_E, Lambda
    elif any(i != None for i in [t,r]):
        aux     = t == None
        if t == None:
            return lambda t1: p_2_E(t1,r1=r), Lambda
        else:
            return lambda r1: p_2_E(t1=t,r1=r1), Lambda
    else:
        return p_2_E(t,r), Lambda

def importData(simulation:str, probe:int=2, time:float = None)-> tuple:
    
    assert 0<= probe <= 10, 'ERRO: O valor da probe deve estar entre 0 e 10'
    PATH_DATA   = Path().absolute().parent / 'data'
    PROBES      = PATH_DATA / 'probes' / simulation /'p.txt'
    FWH         = PATH_DATA / 'acousticData' / simulation / 'FWH-time.dat'
    FWH2        = PATH_DATA / 'acousticData' / simulation / 'FWH-time.dat'

    toPa        = 101325

    if time == None:
        t, p            = np.loadtxt(PROBES, usecols = (0, probe+1), unpack=True)
        fwh_t, fwh_p    = np.loadtxt(FWH, usecols=(0,probe+1), skiprows=1, unpack=True)
        fwh2_t, fwh2_p  = np.loadtxt(FWH2, usecols=(0,probe+1), skiprows=1, unpack=True)
        print(f'Probe: {probe}')
        return ((t,p-toPa), (fwh_t, fwh_p), (fwh2_t,fwh2_p)) 
    else:
        tp      = np.loadtxt(PROBES, usecols = 0)
        pos     = np.searchsorted(tp,time)-1
        p       = np.loadtxt(PROBES, skiprows=12 + pos)[1,1:]

        tfwh    = np.loadtxt(FWH, usecols=0, skiprows=1)
        pos     = np.searchsorted(tfwh,time) - 1
        pfwh    = np.loadtxt(FWH,  skiprows=1+pos)[1:]
        pfwh2   = np.loadtxt(FWH2, skiprows=1+pos)[1:]

        print(f'Time = {tp[pos]} \nPos = {pos}')
        return (p - toPa, pfwh, pfwh)

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
    pfunc,_ = pressure(r = robs)
    aux1    = (efectivePressure - pfunc(efectiveTime))**2
    num     = trapz(aux1, efectiveTime)
    den     = trapz(pfunc(efectiveTime)**2, efectiveTime)
    rms     = num/den

    return rms  

def rmsSpacial(rp:tuple, tobs:float = 0.5, plot:bool = True):
    
    r, p    = rp
    pfunc,_ = pressure(t = tobs)
    
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

def plotTime(FWH:tuple, FWH2:tuple, SIM:tuple, robs:float = None, title:str = None) -> None:
    fwh_t, fwh_p        = FWH
    fwh2_t, fwh2_p      = FWH2
    t, p                = SIM
    
    if robs != None:
        pfunc,_         = pressure(r = robs) 
        plt.plot(t, pfunc(t), 'k', label = 'Analítico', alpha = 0.5)
        plt.plot(t, p,'r-.', label = 'Cálculo Direto')
    else:
        plt.plot(t, p,'r-.', label = 'FWH2')
    plt.plot(fwh_t, fwh_p, 'b--', label = 'FWH', alpha = 0.75)
    plt.plot(fwh2_t, fwh2_p, 'g--', label = 'FWH2', alpha = 0.75)

    if title!=None: plt.title(title)
    plt.xlabel('Tempo [s]')
    plt.ylabel('Pressão [Pa]')
    plt.legend()
    
    plt.grid()
    plt.show()

def plotSpacial(FWH:tuple, SIM:np.ndarray, r:np.ndarray,     
                tobs:float=None, title:str = None)->None:
    pfwh, pfwh2    = FWH
    p              = SIM

    if tobs != None:
        pfunc,_         = pressure(t = tobs)
        rfunc           = np.linspace(r[0],r[-1])
        plt.plot(rfunc, pfunc(rfunc), 'k', label = 'Analítico', alpha =0.5)
        plt.plot(r, p,'ro-.', label = 'Cálculo Direto', alpha = 0.35)
    else:
        plt.plot(r, p,'r-.', label = 'FWH2')
    plt.plot(r, pfwh, 'bo--', label = 'FWH', alpha = 0.75)
    plt.plot(r, pfwh2, 'go--', label = 'FWH', alpha = 0.75)


    if title!=None: plt.title(title)
    plt.xlabel(r'r [m]')
    plt.ylabel('Pressão [Pa]')
    plt.legend()
    
    plt.grid()
    plt.show()



