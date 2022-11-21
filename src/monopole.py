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

from datadrive import PATH_DATA
from pathlib import Path


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
t_num = [2,4]
epsilon_num = p0 * gamma
y_num = 0
xlim = 200
ylim = 200
nx = 401
ny = 401


xrange = np.linspace(-xlim, xlim, nx)
# yrange = np.linspace(-ylim, ylim, ny)
deltax = xrange[1] - xrange[0]
# deltay = yrange[1] - yrange[0]

f = epsilon * sy.exp(-alpha * (x**2 + y**2))

dGdt = sy.diff(G, t)
dGdx = sy.diff(G, x)

H = sy.im(dGdt + M_num * c0 * dGdx)

for j in range(len(t_num)):
    dGdt_num = np.zeros((len(xrange)), dtype=complex)
    dGdx_num = np.zeros((len(xrange)), dtype=complex)
    H_num = np.zeros((len(xrange)), dtype=float)
    f_num = np.zeros((len(xrange)), dtype=float)

    for i in range(len(xrange)):
        H_num[i] = H.xreplace(
            {
                x: xrange[i],
                y: y_num,
                omega: omega_num,
                t: t_num[j],
                c0: c0_num,
                M: M_num,
            }
        ).evalf()
        f_num[i] = f.xreplace(
            {
                x: xrange[i],
                y: y_num,
                alpha: alpha_num,
                epsilon: epsilon_num,
            }
        ).evalf()


    # pprime=-np.imag(sp.signal.fftconvolve(f_num,dGdt_num+M_num*c0_num*dGdx_num,'same')*deltax*deltay)
    pprime = sp.signal.fftconvolve(f_num, H_num, 'same') * deltax 
    xran   = np.arange(len(pprime)) * deltax - xlim

    plt.plot(xran, pprime)
    plt.xlim([-100, 100])
    plt.show()

    arqsave = np.column_stack((xran, pprime))
    pathsave = Path(Path().parent, PATH_DATA, f'teste1{t_num[j]}s.dat').absolute()

    np.savetxt(pathsave, arqsave, header='x P')

# [X, Y] = np.meshgrid(xrange, yrange)
# fig, ax = plt.subplots(1, 1)
# ax.contourf(X, Y, pprime.T)
# ax.set_title('Contour Plot')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_xlim([-100, 100])
# ax.set_ylim([-100, 100])

# plt.show()

# %%
