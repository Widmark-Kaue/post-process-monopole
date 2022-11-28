#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 19:07:26 2022

@author: jsalazar
"""
#%%
import numpy as np
import scipy as sp
import sympy as sy
import matplotlib.pyplot as plt
# from spyder_kernels.utils.iofuncs import save_dictionary


def createDict(*args):
    return dict(((k, eval(k)) for k in args))


# defining variables
x, y, omega, t, c0, M, alpha, epsilon = sy.symbols(
    'x y omega t c0 M alpha epsilon', real=True
)

G = (
    sy.I
    / (4 * c0**2 * sy.sqrt(1 - M**2))
    * sy.hankel1(
        0,
        omega * sy.sqrt(x**2 + (1 - M**2) * y**2) / (c0 * (1 - M**2)),
    )
    * sy.exp(-sy.I * M / (1 - M**2) * omega / c0 * x - sy.I * omega * t)
)

gamma = 1.4
T0 = 298.15
p0 = 101325
R = 8314.46261815324 / 28.9   # Air
omega_num = 20 * np.pi
c0_num = np.sqrt(gamma * R * T0)
alpha_num = np.log(2) / 9
M_num = 0.5
t_num = 2
epsilon_num = p0 * gamma
y_num = 0
xlim = 200
ylim = 5
nx = 401
ny = 11

f = epsilon * sy.exp(-alpha * (x**2 + y**2))

dGdt = sy.diff(G, t)
dGdx = sy.diff(G, x)


xrange = np.linspace(-xlim, xlim, nx)
yrange = np.linspace(-ylim, ylim, ny)
deltax = xrange[1] - xrange[0]
deltay = yrange[1] - yrange[0]

dGdt_num = np.zeros((len(xrange), len(yrange)), dtype=complex)
dGdx_num = np.zeros((len(xrange), len(yrange)), dtype=complex)
f_num = np.zeros((len(xrange), len(yrange)), dtype=float)
H_num = np.zeros((len(xrange), len(yrange)), dtype=float)


H = sy.im(dGdt + M_num * c0 * dGdx)

for i in range(len(xrange)):
    for j in range(len(yrange)):
        H_num[i, j] = H.xreplace(
            {
                x: xrange[i],==0
                y: yrange[j],
                omega: omega_num,
                t: t_num,
                c0: c0_num,
                M: M_num,
            }
        ).evalf()
        f_num[i, j] = f.xreplace(
            {
                x: xrange[i],
                y: yrange[j],
                alpha: alpha_num,
                epsilon: epsilon_num,
            }
        ).evalf()


# pprime=-np.imag(sp.signal.fftconvolve(f_num,dGdt_num+M_num*c0_num*dGdx_num,'same')*deltax*deltay)
pprime = sp.signal.fftconvolve(f_num, H_num, 'same') * deltax * deltay

np.savetxt('pprime.dat', pprime, header='X      Y')

plt.plot(np.arange(len(pprime)) * deltax - xlim, pprime[:, int((ny - 1) / 2)])
plt.xlim([-100, 100])

[X, Y] = np.meshgrid(xrange, yrange)
fig, ax = plt.subplots(1, 1)
ax.contourf(X, Y, pprime.T)
ax.set_title('Contour Plot')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])

plt.show()

data = createDict(
    'X',
    'Y',
    'pprime',
    'gamma',
    'T0',
    'p0',
    'R',
    'omega_num',
    'c0_num',
    'alpha_num',
    'M_num',
    't_num',
    'epsilon_num',
    'y_num',
    'xlim',
    'ylim',
    'nx',
    'ny',
)
outputfile = 'monopole' + '_' + str(M_num) + '.spydata'
# save_dictionary(data, outputfile)
# %%
