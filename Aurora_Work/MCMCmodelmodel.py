import numpy as np

import matplotlib.pyplot as plt

from itertools import chain

from scipy.fft import fft, fftfreq
from scipy import signal
from scipy import constants
import pytensor.tensor as pt
import pymc as pm
mu0 = constants.mu_0
from scipy.special import gamma, jv
import copy

#powerspec,omega, conducitivies,Alfven_speed = np.genfromtxt('powerspectrum.csv', delimiter=','), np.genfromtxt('frequencies.csv', delimiter=',' ), np.genfromtxt('conducitivies.csv', delimiter=','), np.genfromtxt('Alfven.csv', delimiter=',')
#print(conducitivies,Alfven_speed)
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
    def __init__(self, SigmaP, Vai, z, h, epsilon, interval):
        self.SigmaP = SigmaP
        self.Vai = Vai
        self.z = z
        self.h = h
        self.epsilon = epsilon
        self.SigmaA = 1 / (
            Vai * mu0
        )  # Alven conductance as given by equation on the right side at the top of pg 1556
        self.omega  = omega * 2 * np.pi

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


def my_loglike(theta, data):
    Vai=theta[1]*5e5
    model_init = model(
        theta[0],
        Vai,
        theta[3], #Z
        theta[2]*120e3, #H
        theta[4],
        0.01
    )

    dat = model_init.r_bayes()[1]
    
    model_data = np.abs(dat * Vai)

    return -1/(2*theta[5]**2) * (np.log10(data) - np.log10(model_data)) ** 2 - np.log(np.sqrt(2 * np.pi)) - np.log(theta[5]) # Simple distance minimization


# define a pytensor Op for our likelihood function
class LogLike(pt.Op):

    """
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """

    itypes = [pt.dvector]  # expects a vector of parameter values when called
    otypes = [pt.dvector]  # outputs a single scalar value (the log likelihood)

    def __init__(self, loglike, data):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that our function requires.
        """

        # add inputs as class attributes
        self.likelihood = loglike
        self.data = data

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        (theta,) = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta, self.data)

        outputs[0][0] = np.array(logl)  # output the log-likelihood

#TODO need height
omega = np.array([
    0.080, 0.160, 0.240, 0.320, 0.400, 0.480, 0.560, 0.640, 0.720, 0.800,
    0.880, 0.960, 1.04, 1.12, 1.20, 1.28, 1.36, 1.44, 1.52, 1.60, 1.68, 1.76,
    1.84, 1.92, 2.00, 2.08, 2.16, 2.24, 2.32, 2.40, 2.48, 2.56, 2.64, 2.72,
    2.80, 2.88, 2.96, 3.04, 3.12, 3.20, 3.28, 3.36, 3.44, 3.52, 3.60, 3.68,
    3.76, 3.84, 3.92, 4.00, 4.08, 4.16, 4.24, 4.32, 4.40, 4.48, 4.56, 4.64,
    4.72, 4.80, 4.88, 4.96, 5.04, 5.12, 5.20, 5.28, 5.36, 5.44, 5.52, 5.60,
    5.68, 5.76, 5.84, 5.92, 6.00
])
model_init = model(
    1.35,
    3.532e5,
    328e3,
    300e3,
    0.01,
    0.01)

dat = model_init.r_bayes()[1]
powerspec = np.abs(dat * 3.532e5)
conducitivies=2
def main():
    logl = LogLike(my_loglike, powerspec)
    with pm.Model():
        # uniform priors on m and c
        
        sigmap = pm.Uniform("sigmap", lower=0.1, upper=12, initval=conducitivies)

        VA_norm = pm.Uniform("va", lower=0.3, upper=3, initval=1)
        h_norm = pm.Uniform("h", lower=0.3, upper=3, initval=1)
        deviation= pm.HalfCauchy("deviation", beta=10 ) # from https://www.pymc.io/projects/docs/en/stable/learn/core_notebooks/GLM_linear.html#glm-linear
        # convert m and c to a tensor vector
        theta = pt.as_tensor_variable([sigmap, VA_norm, h_norm, 428e3-100e3, 0.01, deviation])

        # use a Potential to "call" the Op and include it in the logp computation
        pm.Potential("likelihood", logl(theta))

        # Use custom number of draws to replace the HMC based defaults
        idata_mh = pm.sample(5, tune=5, cores=6, return_inferencedata=True)
    idata_mh.to_netcdf("convergencetest1.nc")


if __name__ == "__main__":
    main()


