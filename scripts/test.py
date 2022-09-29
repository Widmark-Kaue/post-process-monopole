#%% Librarys
import numpy as np
import wave as wv

#%% Parâmetros de simulação
microphone = np.arange(2, 112, 10)

#%% Comportamento temporal
# # SIM, FWH40 = importData('backward40')
# _,FWH20    = importData('backward20')

for m in range(len(microphone)):
    SIM,FWH, FWH2      = wv.importData('monopole3', probe = m)


    #% Plot
    ti = f'Comportamento temporal (r = {microphone[m]} [m])'
    wv.plotTime(FWH, FWH2 ,SIM, robs=microphone[m], title = ti)

#%% Comportamento espacial
p, pfwh, pfwh2 = wv.importData('monopole3', time  = 0.5)
ti = f'Comportamento espacial (t = 0.5 [s])'
wv.plotSpacial((pfwh,pfwh2), p, r = microphone, tobs=0.5, title=ti)

#%% RMS
RMS = wv.rmsTime(SIM, r = microphone[2])*100
print(RMS)

