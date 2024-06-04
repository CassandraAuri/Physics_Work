import numpy as np
import matplotlib.pyplot as plt
from PyEMD import EMD

# Generate a sample non-stationary signal
np.random.seed(0)
t = np.linspace(0, 1, 1000)
signal = np.sin(2 * np.pi * 5 * t**2) + np.sin(2 * np.pi * 30 * t)

# Perform EMD
emd = EMD()
imfs = emd(signal)

# Plot IMFs
plt.figure()
for i, imf in enumerate(imfs):
    plt.subplot(len(imfs), 1, i + 1)
    plt.plot(t, imf)
    plt.title(f'IMF {i + 1}')
plt.tight_layout()
plt.show()