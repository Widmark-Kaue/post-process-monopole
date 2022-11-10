#%% Librarys
import numpy as np
import matplotlib.pyplot as plt

from re import findall
from pathlib import Path
from scipy.integrate import trapz
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
    t: float or np.ndarray,
    y: float or np.ndarray = 0,
    M: float = 0.1,
):

    assert 0 <= M <= 1, 'Número de Mach deve está entre 0 e 1'
    if type(y) != np.ndarray:
        y = y * np.ones(len(x))

    # Constantes Globais
    c0 = 340.29
    c = 331.45
    T0 = 273.15
    T = (c0 / c) ** 2 * T0
    rho0 = 101325 / (287.058 * T)   #% Eq. dos Gases ideias
    S = 0.1
    freq = 100
    omega = freq * 2 * np.pi
    k = omega / c0

    # Termos da solução (eqs. 4.11 4.13 )
    csi = (
        omega
        * np.sqrt(x**2 + (1 - M**2) * y**2)
        / ((1 - M**2) * c0**2)
    )

    eta = -1j * M * k * x / (1 - M**2) - 1j * omega * t
    H0Flow = hankel1(0, csi)
    H1Flow = hankel1(1, csi)

    T1 = omega / (4 * c0**3 * (1 - M**2) ** (3 / 2))
    T2 = M * H0Flow - 1j * x * H1Flow / np.sqrt(x**2 + (1 - M**2) * y**2)
    T3 = np.exp(eta)
    Gx = T1 * T2 * T3   # derivada em x da função de Green
    Gt = (
        1j * (-1j * omega) / (4 * c0**2 * np.sqrt(1 - M**2)) * H0Flow * T3
    )   # derivada em t da função de Green
    f = rho0 * (c0**2) * S
    pFlow = (f * (Gt + M * Gx)).imag

    return pFlow

def pressureFlow2(
    xlim:tuple  = (200,-200),
    ylim:tuple  = (200,-200),
    nxy:tuple   = (401,401),
    t: float    = 0.6,
    alpha:float = np.log(2)/9,
    freq:float  = 10,
    M: float    = 0.5, 
    gamma:float = 1.4,
    PTR: tuple  = (101325, 298.15, 8314.46261815324/28.9 )   
):

    assert 0 <= M <= 1, 'Número de Mach deve está entre 0 e 1'

    # vetores
    nx, ny  = nxy 
    x       = np.linspace(xlim[1], xlim[0], nx)
    y       = np.linspace(ylim[1],ylim[0], ny)
    deltax  = x[1] - x[0]
    deltay  = y[1] - y[0]

    f       = np.zeros(nxy, dtype=float)
    H       = np.zeros(nxy, dtype=float)
    

    # Constantes Globais
    p0, T0, R   = PTR    
    omega       = freq * 2 * np.pi
    c0          = np.sqrt(gamma * R * T0)
    k           = omega/c0
    epsilon     = p0 * gamma

    for i,xi in enumerate(x):
        for j,yj in enumerate(y):
            try:
                del(G, dGdt, dGdx)
            except:
                pass

            # Gaussian distribution of amplitude
            f[i,j] = epsilon*np.exp(-alpha*(xi**2+yj**2))

            # Green's function
            if xi == 0: 
                xksi = 1e-100 # evitar indeterminação
            else:
                xksi = xi
            ksi = omega*np.sqrt(xksi**2 + (1 - M**2) * yj**2)/((1 - M**2) * c0)
            eta = -1j*M*k*xi/(1 - M**2)-1j*omega*t

            G = 1j/(4*c0**2 + np.sqrt(1 - M**2))
            G*= hankel1(0,ksi)*np.exp(eta)

            # Derivate of Green's function
            dGdx = omega/(4*c0**3 * (1-M**2)**(3/2))
            if xi == 0:
                dGdx*=(M*hankel1(0, ksi)-1j*hankel1(1,ksi))
            else:
                dGdx*=(M*hankel1(0, ksi)-1j*xi*hankel1(1,ksi)/np.sqrt(xi**2+(1-M**2)*yj**2))
            dGdx*= np.exp(eta)

            dGdt = G*(-1j*omega)

            # Solution the convolution product
            H[i,j] = (dGdt + M*dGdx).imag
    
    pFlow = fftconvolve(f, H, 'same') * deltax * deltay
    print(len(pFlow))
    print(pFlow.shape)

    # Plot
    print(ny)
    Xplt = np.arange(len(pFlow))*deltax - abs(xlim[1])
    Yplt = pFlow[:,int( (ny -1)/2 )] 
    plt.plot(Xplt, Yplt)
    plt.xlim([-100,100])
    plt.xlabel('x [m]')
    plt.ylabel('P [Pa]')
    plt.grid()
    plt.show()

    [X,Y] = np.meshgrid(x,y)
    fig, ax = plt.subplots(1,1)
    ax.contourf(X, Y, pFlow.T)
    ax.set_title('Contour Plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-100, 100])
    ax.set_ylim([-100, 100])

    plt.show()

    return pFlow


def importData(
    simulation: str, probe: int = 2, time: float = 0, case: str = 'monopole'
) -> tuple:

    PATH_DATA = Path(Path().absolute().parent, 'data', case)
    PROBES = Path(PATH_DATA, 'probes', simulation, str(time), 'p.txt')
    FWH = Path(PATH_DATA, 'acousticData', simulation, 'FWH-time.dat')
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
