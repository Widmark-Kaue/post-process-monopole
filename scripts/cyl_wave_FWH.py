#%% Librarys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel2
from pathlib import Path
#%% Caminhos
PATH_DATA   = Path('/home/labcc/Github/post-process-monopole/data')
PROBES      = PATH_DATA / 'probes' / '0'
FWH         = PATH_DATA / 'pressureData' 

#%% Constantes
rvec = np.arange(2, 102, 10)

raio = 0.05715 / 2
S = 0.1
c0 = 340.29
c = 331.45
freq = 100
omega = freq * 2 * np.pi
T0 = 273.15
T = (c0 / c) ** 2 * T0
#%% Cálculo da densidade
rho0 = 101325 / (287.058 * T)
area = 2 * np.pi * raio
velocity = S / area
H_1_fonte_J = hankel2(1, (omega * raio / c0))

A0 = velocity*1j*rho0*c0/H_1_fonte_J

#%% Import dados
probe        = 2
t, p         = np.loadtxt(PROBES / 'p.txt', 
                          usecols = (0, probe-1), unpack=True)
fwh_t, fwh_p = np.loadtxt(FWH / 'FWH_surfacePressure.txt', 
                           usecols=(0,probe-1), unpack=True)
 

#%% Analítico
H_0_E = hankel2(0, omega*rvec/c0)
P_2_E = A0*H_0_E
p_2_E = P_2_E*np.exp(1j*omega*t).imag

# %% Plot
plt.plot(t, p_2_E[0:len(t)], 'k--', label = 'Analítico')
plt.plot(t, p - 101325,'r--', label = 'Cálculo Direto')
plt.plot(fwh_t, fwh_p, 'b--', label = 'FWH')

plt.xlabel('Tempo [s]')
plt.ylabel('Pressão [Pa]')
plt.legend()

plt.show()
# %%
