#%% Librarys
import numpy as np
import wave as wv
#%% Parâmetros de simulação
microphone = np.linspace(2, 102, 20)

#%% Comportamento temporal
for m in range(len(microphone)):
    SIM, FWH, FWH2 = wv.importData('32ppw', probe=m)

    #% Plot
    ti = f'Comportamento temporal (r = {round(microphone[m],2)} [m])'
    wv.plotTime(FWH, FWH2, SIM, robs=microphone[m], title=ti)

#%% Comportamento espacial
p, pfwh, pfwh2 = wv.importData('back40', time=0.5)
ti = f'Comportamento espacial (t = 0.5 s)'
wv.plotSpacial((pfwh, pfwh2), p, r=microphone, tobs=0.5, title=ti)

#%% RMS
RMS = wv.rmsTime(SIM, robs=microphone[0]) * 100
print(RMS)


# %% teste



# %%
