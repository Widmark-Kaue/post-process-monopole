#%% Librarys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hankel1, hankel2
from pathlib import Path

PATH_DATA = Path('/home/labcc/Github/post-process-monopole/data')
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
rho = 101325 / (287.058 * T)
area = 2 * np.pi * raio
velocity = S / area
H_1_fonte_J = hankel2(1, (omega * raio / c0))

A0 = velocity*1j*rho0*c0/H_1_fonte_J

#%% Import dados


#%% Analítico
