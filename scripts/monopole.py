#%% Librarys
import numpy as np
import wave as wv

#%% Parâmetros de simulação
microphone = np.linspace(2, 102, 11)
"""
32 PPW
"""
"""  
T1
"""
#%% Comportamento temporal
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('32ppwT1', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m])'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)

#%% Comportamento espacial
p, pfwh, pfwh2 = wv.importData('32ppwT1', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s) - n = 4000'
wv.plotSpacial((pfwh, pfwh2), p, r=(2, 102), tobs=0.5, title=ti)

#%% RMS
RMS = wv.rmsTime(SIM, robs=microphone[0]) * 100
print(RMS)
#%%
"""
T7 - n = 80
"""
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('32ppwT7', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m]) - n = 80'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)


p, pfwh, pfwh2 = wv.importData('32ppwT7', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s) - n = 80'
wv.plotSpacial((pfwh, pfwh2), p, r=(2, 102), tobs=0.5, title=ti)

"""
T8 - n = 40
"""
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('32ppwT8', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m]) - n = 40'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)

p, pfwh, pfwh2 = wv.importData('32ppwT8', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s) - n = 40'
wv.plotSpacial((pfwh, pfwh2), p, r=(2, 102), tobs=0.5, title=ti)
#%%
""" 
64 PPW
"""
"""
T8 - n = 80
"""
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('64ppwT8', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m]) - n = 80'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)


p, pfwh, pfwh2 = wv.importData('64ppwT8', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s) - n = 80'
wv.plotSpacial((pfwh, pfwh2), p, r=(2, 102), tobs=0.5, title=ti)

"""
T9 - n = 40
"""
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('64ppwT9', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m]) - n = 40'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)

p, pfwh, pfwh2 = wv.importData('64ppwT9', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s) - n = 40'
wv.plotSpacial((pfwh, pfwh2), p, r=(2, 102), tobs=0.5, title=ti)
# %%
