import numpy as np
def getspectrum(min_mz, max_mz, n_peaks):
    return min_mz + max_mz*np.random.rand(n_peaks), np.abs(np.random.randn(n_peaks))