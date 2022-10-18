#%% Library
import numpy as np
import wave as wv
import matplotlib.pyplot as plt

#%% Teste 3
pSim,pFWH,pFWH2 = wv.importData('pimple', time = 0.5, case = 'monopoleFlow')
x = np.linspace(2,102,len(pSim))

plt.plot(x,pSim, label = 'Cálculo Direto')
plt.xlabel('x [m]')
plt.ylabel('p [Pa]')
plt.title('Teste ')
plt.grid()
plt.legend()
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
