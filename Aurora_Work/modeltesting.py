# inttatlize packages
import numpy as np

# import mplcyberpunk

from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from scipy import constants
mu0 = constants.mu_0
from scipy.special import gamma, jv
import matplotlib.animation as animation
import copy
import numpy as np
import pandas as pd
filename = 'powerspectrum.csv'
data = pd.read_csv(filename, sep=",", header=None)

def cbesselj2(nu, x):
    kMax = 100
    tol = 1e-2
    x = np.array(x).reshape(-1, 1)
    Jnu = np.zeros_like(x, dtype=complex)
    for m in range(len(x)):
        if np.isreal(nu[m]):
            Jnu[m] = jv(nu[m], x[m])
        else:
            if x[m] < 30:
                Jnu[m] = 0.0 + 0j
                for k in range(kMax + 1):
                    JnuLast = copy.deepcopy(Jnu[m])
                    Jnu[m] += (-0.25 * x[m] * x[m]) ** k / (
                        gamma(k + 1) * gamma(nu[m] + k + 1)
                    )
                    if np.real(Jnu[m]) != 0:
                        Rerr = abs(
                            (np.real(Jnu[m]) - np.real(JnuLast)) / np.real(Jnu[m])
                        )
                    if np.imag(Jnu[m]) != 0:
                        Ierr = abs(
                            (np.imag(Jnu[m]) - np.imag(JnuLast)) / np.imag(Jnu[m])
                        )
                    ##print((Rerr), (Ierr))
                    if "Ierr" in locals() and "Rerr" in locals():
                        if np.max(Rerr) <= tol and np.max(Ierr) <= tol:
                            ##print("real")
                            break
                    elif "Rerr" in locals():
                        if np.max(Rerr) <= tol:
                            ##print("false")
                            break
                Jnu[m] = Jnu[m] * (0.5 * x[m]) ** nu[m]
                if np.real(k) == kMax:
                    print(
                        "Algorithm does not converge in the calculation of bessel function. Maximum concurence number arrived!"
                    )
            else:
                Jnu[m] = np.sqrt(2 / (np.pi * x[m])) * np.cos(
                    x[m] - nu[m] * (np.pi / 2) - (np.pi / 4)
                )
    Jnu = Jnu.flatten()
    return Jnu


class model(object):
    def __init__(self, SigmaP, Vai, h, epsilon, interval):
        self.SigmaP = SigmaP
        self.Vai = Vai
        self.z = 330e3
        self.h = h
        self.epsilon = epsilon
        self.SigmaA = 1 / (
            Vai * mu0
        )  # Alven conductance as given by equation on the right side at the top of pg 1556
        self.omega = np.arange(interval, 8 + interval, interval) * 2 * np.pi

    def x0(self):
        return 2 * self.h * self.omega / self.Vai

    def x(self):
        return self.x0() * np.exp(-self.z / (2 * self.h))

    def r_bayes(self):
        alpha = mu0 * self.SigmaP * self.Vai
        x0 = self.x0()
        x = self.x()
        epsilon = self.epsilon

        intervalx0 = x0[1] - x0[0]
        intervalx = x[1] - x[0]

        ##R = np.zeros(np.shape(x0)[0], dtype=complex)
        ##phiOverVaiAz = np.zeros(np.shape(x0)[0], dtype=complex)
        ##normalisedPhiOverVaiAz = np.zeros(np.shape(x0)[0]).astype(complex)
        ##phiOverAz = np.zeros(np.shape(x0)[0]).astype(complex)
        ##impedance = np.zeros(np.shape(x0)[0]).astype(complex)

        num2 = cbesselj2(1j * x0 * epsilon, x0)
        denom2 = cbesselj2(-1j * x0 * epsilon, x0)

        rationum1 = cbesselj2(1j * x0 * epsilon, x)
        rationum2 = cbesselj2(-1j * x0 * epsilon, x)

        num2tmp1 = cbesselj2(1j * (x0) * epsilon, (x0 - (intervalx0 / 2)))
        denom2tmp1 = cbesselj2(-1j * (x0) * epsilon, (x0 - (intervalx0 / 2)))
        rationum1tmp1 = cbesselj2(1j * (x0) * epsilon, (x - (intervalx / 2)))
        rationum2tmp1 = cbesselj2(-1j * (x0) * epsilon, (x - (intervalx / 2)))

        num2tmp2 = cbesselj2(1j * (x0) * epsilon, (x0 + (intervalx0 / 2)))
        denom2tmp2 = cbesselj2(-1j * (x0) * epsilon, (x0 + (intervalx0 / 2)))
        rationum1tmp2 = cbesselj2(1j * (x0) * epsilon, (x + (intervalx / 2)))
        rationum2tmp2 = cbesselj2(-1j * (x0) * epsilon, (x + (intervalx / 2)))

        num1 = (num2tmp2 - num2tmp1) / intervalx0
        denom1 = (denom2tmp2 - denom2tmp1) / intervalx0
        ratiodenom1 = (rationum1tmp2 - rationum1tmp1) / intervalx
        ratiodenom2 = (rationum2tmp2 - rationum2tmp1) / intervalx

        R = -(1j * num1 + alpha * num2) / (1j * denom1 + alpha * denom2)
        phiOverVaiAz = (
            -1j
            * (x0 / x)
            * ((rationum1 + R * rationum2) / (ratiodenom1 + R * ratiodenom2))
        )
        # normalisedPhiOverVaiAz[i] = phiOverVaiAz*alpha
        # phiOverAz[i] = phiOverVaiAz*self.Vai
        # impedance[i] = phiOverAz*mu0
        return R, np.abs(phiOverVaiAz)



model_init = model(
    1,
    1e6*0.2,
    100e3*1.5,
    0.01,
    0.15)

dat = model_init.r_bayes()[1]
model_data = np.abs(dat * 1e6)
plt.plot(np.log10(model_data), color='orange')
plt.plot(np.log10(data))
plt.show()

