#%% Librarys
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

from re import findall
from pathlib import Path
from scipy.integrate import trapz
from src.datadrive import PATH_DATA
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
from scipy.special import hankel2, hankel1


def pressure(
    r: float or np.ndarray = None,
    t: float or np.ndarray = None,
) -> float or np.ndarray or any:

    # constantes globais
    rf = 0.05715 / 2
    S = 0.1
    c0 = 340.29
    c = 331.45
    freq = 100
    omega = freq * 2 * np.pi
    T0 = 273.15
    T = (c0 / c) ** 2 * T0
    rho0 = 101325 / (287.058 * T)   #% Eq. dos Gases ideias
    area = 2 * np.pi * rf
    velocity = S / area
    Lambda = c / freq

    H_1_fonte_J = hankel2(1, (omega * rf / c0))
    A0 = velocity * 1j * rho0 * c0 / H_1_fonte_J

    H_0_E = lambda r1: hankel2(0, omega * r1 / c0)
    P_2_E = lambda r1: A0 * H_0_E(r1)
    p_2_E = lambda t1, r1: (P_2_E(r1) * np.exp(1j * omega * t1)).imag

    if r == None and t == None:
        return p_2_E, Lambda
    elif any(i != None for i in [t, r]):
        aux = t == None
        if t == None:
            return lambda t1: p_2_E(t1, r1=r), Lambda
        else:
            return lambda r1: p_2_E(t1=t, r1=r1), Lambda
    else:
        return p_2_E(t, r), Lambda


def pressureFlow(
    x: float or np.ndarray,
    y: float or np.ndarray = 0,
    t: float = 0.6,
    alpha: float = np.log(2) / 9,
    freq: float = 10,
    M: float = 0.5,
    gamma: float = 1.4,
    PTR: tuple = (101325, 298.15, 8314.46261815324 / 28.9),
):

    assert 0 <= M <= 1, 'Número de Mach deve está entre 0 e 1'

    # Constantes Globais
    p0, T0, R = PTR
    omega = freq * 2 * np.pi
    c0 = np.sqrt(gamma * R * T0)
    k = omega / c0
    epsilon = p0 * gamma
    deltax = x[1] - x[0]

    # Gaussian distribution of amplitude
    f = epsilon * np.exp(-alpha * (x**2 + y**2))

    # Green's function
    ksi = omega * np.sqrt(x**2 + (1 - M**2) * y**2) / ((1 - M**2) * c0)
    eta = -1j * M * k * x / (1 - M**2) - 1j * omega * t

    G = 1j / (4 * c0**2 + np.sqrt(1 - M**2))
    G *= hankel1(0, ksi) * np.exp(eta)

    # Derivate of Green's function
    dGdx = omega / (4 * c0**3 * (1 - M**2) ** (3 / 2))
    if y == 0 and (0 in x):
        x[list(x).index(0)] = 1

    dGdx *= M * hankel1(0, ksi) - 1j * x * hankel1(1, ksi) / np.sqrt(
        x**2 + (1 - M**2) * y**2
    )
    dGdx *= np.exp(eta)

    dGdt = G * (-1j * omega)

    # Solution for the convolution product
    H = np.imag(dGdt + M * c0 * dGdx)

    pFlow = fftconvolve(f, H, 'same') * deltax

    return pFlow


def pressureFlow2(
    xlim: tuple = (-200, 200),
    ylim: tuple = (-200, 200),
    nxy: tuple = (401, 401),
    t: float = 0.6,
    alpha: float = np.log(2) / 9,
    freq: float = 10,
    M: float = 0.5,
    gamma: float = 1.4,
    PTR: tuple = (101325, 298.15, 8314.46261815324 / 28.9),
):

    assert 0 <= M <= 1, 'Número de Mach deve está entre 0 e 1'

    # Vetores
    nx, ny  = nxy
    x       = np.linspace(xlim[0], xlim[1], nx, dtype=float)
    y       = np.linspace(ylim[0], ylim[1], ny, dtype=float)
    deltax  = x[1] - x[0]
    deltay  = y[1] - y[0]

    # Funções da convolução
    f = np.zeros(nxy, dtype=float)
    H = np.zeros(nxy, dtype=float)


    # refino em x
    
    # search_pos = lambda arr, val: int(np.argwhere(arr==val)[0][0])
    # if 0 in x:
    #     n0    = search_pos(x,0)
    #     x     = np.delete(x, n0)
    #     y     = np.delete(y, n0)
        # xIns  = np.linspace(x[n0-1],x[n0+1],int(0.2*nx))
        # xIns  = np.delete(xIns, (0,-1))
        # if 0 in xIns:
        #     n1 = search_pos(xIns, 0)
        #     xIns = np.delete(xIns, n1)
        
        # x = np.delete(np.insert(x,n0, xIns), n0+len(xIns))
        # print(f'tam x = {len(x)}')
        # print(f'tam y = {len(y)}')
        # assert not (0 in x), "puta que pariu viu!"

   

    # Constantes Globais
    p0, T0, R   = PTR
    omega       = freq * 2 * np.pi
    c0          = np.sqrt(gamma * R * T0)
    k           = omega / c0
    epsilon     = p0 * gamma

    # Variáveis símbolicas
    xSy, ySy, tSy, etaSy = sy.symbols('xSy ySy tSy etaSy', real = True)

    GSy = (
        sy.I/ (4 * c0**2 * sy.sqrt(1 - M**2))
        * sy.hankel1(
            0,
            omega * sy.sqrt(xSy**2 + (1 - M**2) * ySy**2) 
            / (c0 * (1 - M**2)))
        * sy.exp(etaSy)
        )

    hankel01Sy = (
        sy.hankel1(
                0,
                omega*sy.sqrt(xSy**2+(1-M**2)*ySy**2)
                /(c0*(1- M**2))
                  )
                )
    hankel11Sy = (
        sy.hankel1(
            1,
            omega*sy.sqrt(xSy**2+(1-M**2)*ySy**2)
            /(c0*(1- M**2))
                  )
                )

    dGdxSy = omega / (4 * c0**3 * (1 - M**2) ** (3 / 2))
    dGdxSy *= (
        M * hankel01Sy- sy.I * (xSy) * hankel11Sy
        / sy.sqrt(xSy**2 + (1 - M**2) * ySy**2) 
            )
    dGdxSy *= sy.exp(etaSy)

    dGdtSy = GSy * (-sy.I * omega)
    # dGdtSy = sy.diff(GSy, tSy)
    # dGdxSy = sy.diff(GSy, xSy)

    HSy = sy.im(dGdtSy + M * c0 * dGdxSy)

    for i, xi in enumerate(x):
        if i%100 ==1:
            print(f'H[{i-1}] =\n{H[i-1,:11]}')
        for j, yj in enumerate(y):
            try:
                del (G, dGdt, dGdx)
            except:
                pass

            # Gaussian distribution of amplitude
            f[i, j] = epsilon * np.exp(-alpha * (xi**2 + yj**2))

            # Green's function
            ksi = (omega* np.sqrt(xi**2 + (1 - M**2) * yj**2)/ ((1 - M**2) * c0))
            eta = -1j * M * k * xi / (1 - M**2) - 1j * omega * t

            if xi == 0 and yj:
                H[i,j] = HSy.xreplace(
                    {
                        xSy: xi,
                        ySy: yj,
                        etaSy: eta,
                    }
                ).evalf()
                print('Passou aqui')
            else:
                G = 1j / (4 * c0**2 + np.sqrt(1 - M**2))
                G *= hankel1(0, ksi) * np.exp(eta)

                # Derivate of Green's function
                dGdx = omega / (4 * c0**3 * (1 - M**2) ** (3 / 2))
                dGdx *= (
                    M * hankel1(0, ksi) - 1j * xi * hankel1(1, ksi) 
                    / np.sqrt(xi**2 + (1 - M**2) * yj**2)
                        )
                dGdx *= np.exp(eta)

                dGdt = G * (-1j * omega)

                # Solution the convolution product
                H[i, j] = (dGdt + M * c0 * dGdx).imag

    pFlow = fftconvolve(f, H, 'same') * deltax * deltay
    print(f'pFlow =\n{pFlow[:11, int((ny-1)/2)]}')

    # Plot
    Xplt = np.arange(len(pFlow)) * deltax - abs(xlim[0])
    Yplt = pFlow[:, int((ny - 1) / 2)]
    plt.plot(Xplt, Yplt)
    plt.xlim([-100, 100])
    plt.xlabel('x [m]')
    plt.ylabel('P [Pa]')
    plt.grid()
    plt.show()

    [X, Y] = np.meshgrid(x, y)
    fig, ax = plt.subplots(1, 1)
    ax.contourf(X, Y, pFlow.T)
    ax.set_title('Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-100, 100])
    ax.set_ylim([-100, 100])

    plt.show()

    return pFlow, (Xplt, Yplt)


def importData(
    simulation: str, probe: int = 2, time: float = 0, case: str = 'monopole'
) -> tuple:

    PATH_DATA_CASE = Path(PATH_DATA, case)
    PROBES = Path(PATH_DATA_CASE, 'probes', simulation, str(time), 'p.txt')
    FWH = Path(PATH_DATA_CASE, 'acousticData', simulation, 'FWH-time.dat')
    FWH2 = FWH.with_name('FWH2-time.dat')

    toPa = 101325

    if time == 0:
        t, p = (
            np.loadtxt(PROBES, usecols=(0, probe + 1), unpack=True)
            if PROBES.exists()
            else (0, 0)
        )

        fwh_t, fwh_p = (
            np.loadtxt(FWH, usecols=(0, probe + 1), skiprows=1, unpack=True)
            if FWH.exists()
            else (0, 0)
        )

        fwh2_t, fwh2_p = (
            np.loadtxt(FWH2, usecols=(0, probe + 1), skiprows=1, unpack=True)
            if FWH2.exists()
            else (0, 0)
        )

        print(f'Probe: {probe}')
        return ((t, p - toPa), (fwh_t, fwh_p), (fwh2_t, fwh2_p))
    else:
        if not PROBES.exists():
            PROBES = PROBES.parent.parent / '0' / 'p.txt'

        tsim = np.loadtxt(PROBES, usecols=0, ndmin=1)

        arq = open(PROBES, 'r')
        skip = len(findall('#', arq.read())) - 1
        arq.close()

        pos1 = np.searchsorted(tsim, time)
        p = np.loadtxt(
            PROBES, skiprows=skip + pos1 - 1 * (pos1 != 0), max_rows=1
        )[1:]

        if FWH.exists():
            tfwh = np.loadtxt(FWH, usecols=0, skiprows=1, ndmin=1)
            pos2 = np.searchsorted(tfwh, time)
            pfwh = np.loadtxt(
                FWH, skiprows=1 + pos2 - 1 * (pos2 != 0), max_rows=1
            )[1:]
            pfwh2 = np.loadtxt(
                FWH2, skiprows=1 + pos2 - 1 * (pos2 != 0), max_rows=1
            )[1:]
        else:
            pfwh = pfwh2 = 0
        print(f'Time = {tsim[pos1]} \nPos = {pos1}')
        return (p - toPa, pfwh, pfwh2)


def rmsTime(
    tp: tuple,
    robs: float,
    freq: float = 100,
    c0: float = 340.29,
    plot: bool = True,
) -> float:

    t, p = tp
    transientTime = robs / c0 + 5 / freq
    fil = lambda x: x >= transientTime
    efectiveTime = np.array(list(filter(fil, t)))
    t0 = list(t).index(efectiveTime[0])
    efectivePressure = p[t0:]

    #% RMS
    pfunc, _ = pressure(r=robs)
    aux1 = (efectivePressure - pfunc(efectiveTime)) ** 2
    num = trapz(aux1, efectiveTime)
    den = trapz(pfunc(efectiveTime) ** 2, efectiveTime)
    rms = num / den

    #% Plot
    if plot:
        plt.plot(efectiveTime, pfunc(efectiveTime), 'k-', label='Analítico')
        plt.plot(efectiveTime, efectivePressure, 'r-', label='Direto')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Pressão [Pa]')
        plt.legend()
        plt.grid()
        plt.show()

    return rms


def rmsSpacial(rp: tuple, tobs: float = 0.5, plot: bool = True) -> float:

    r, p = rp
    pfunc, _ = pressure(t=tobs)

    if plot:
        plt.plot(r, pfunc(r), 'k-', label='Analítico')
        plt.plot(r, p, 'r-', label='Direto')
        plt.xlabel('raio [m]')
        plt.ylabel('Pressão [Pa]')
        plt.legend()
        plt.grid()
        plt.show()

    aux1 = (p - pfunc(r)) ** 2
    num = trapz(aux1, r)
    den = trapz(pfunc(r) ** 2, r)
    rms = num / den

    return rms


def phaseAmplitude():
    return 0


def plotTime(
    FWH: tuple, FWH2: tuple, SIM: tuple, robs: float = None, title: str = None
) -> None:
    fwh_t, fwh_p = FWH
    fwh2_t, fwh2_p = FWH2
    t, p = SIM

    if robs != None:
        pfunc, _ = pressure(r=robs)
        plt.plot(t, pfunc(t), 'k', label='Analítico', alpha=0.5)
        plt.plot(t, p, 'r-.', label='Cálculo Direto')

    plt.plot(fwh_t, fwh_p, 'b--', label='FWH', alpha=1)
    plt.plot(fwh2_t, fwh2_p, 'g--', label='FWH2', alpha=0.75)

    if title != None:
        plt.title(title)
    plt.xlabel('Tempo [s]')
    plt.ylabel('Pressão [Pa]')
    plt.legend()

    plt.grid()
    plt.show()


def plotSpacial(
    FWH: tuple,
    SIM: np.ndarray,
    r: tuple = (2, 102),
    tobs: float = 0.5,
    npts: int = 800,
    title: str = None,
) -> None:
    pfwh, pfwh2 = FWH
    p = SIM
    pfunc, Lambda = pressure(t=tobs)
    rfunc = np.linspace(0, r[-1], npts)
    rvec1 = np.linspace(r[0], r[-1], len(p))
    rvec2 = np.linspace(r[0], r[-1], len(pfwh))

    plt.plot(rfunc / Lambda, pfunc(rfunc), 'k', label='Analítico', alpha=0.5)
    plt.plot(rvec1 / Lambda, p, 'r--', label='Direto', alpha=0.85)
    plt.plot(rvec2 / Lambda, pfwh, 'b--', label='FWH', alpha=0.25)
    plt.plot(rvec2 / Lambda, pfwh2, 'g--', label='FWH2', alpha=0.15)

    if title != None:
        plt.title(title)
    plt.xlabel(r'r/$\lambda$ [m]')
    plt.ylabel('Pressão [Pa]')
    plt.ylim([-10, 10])

    plt.legend()
    plt.grid()
    plt.show()

# %%
