#%% Library
import numpy as np
import wave as wv
import matplotlib.pyplot as plt

#%% Teste 3



#%% teste função
x = np.linspace(-100,100,1000)
pFlow = wv.pressureFlow(x = x, t = 0.5)

plt.plot(x,pFlow, label = "analytical")
plt.xlabel("x")
plt.ylabel("p")
plt.grid()
plt.legend()
plt.plot()

# %%
