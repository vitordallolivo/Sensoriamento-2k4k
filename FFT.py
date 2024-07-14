import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd

def fftPlot(sig, dt=None, plot=True):

    t = np.arange(0, len(sig)) * dt
    xLabel = 'freq [Hz]'

    if sig[0] % 2 != 0:
        warnings.warn("signal preferred to be even in size, autoFixing it...")
        t = t[0:-1]
        sig = sig[0:-1]

    sigFFT = np.fft.fft(sig) / t.shape[0]  # Divided by size t for coherent magnitude

    freq = np.fft.fftfreq(t.shape[0], d=dt)

    # Plot analytic signal - right half of frequence axis needed only...
    firstNegInd = np.argmax(freq < 0)
    freqAxisPos = freq[0:firstNegInd]
    sigFFTPos = 2 * sigFFT[0:firstNegInd]  # *2 because of magnitude of analytic signal

    if plot:
        plt.figure()
        plt.plot(freqAxisPos, np.abs(sigFFTPos))
        plt.xlabel(xLabel)
        plt.ylabel('mag')
        plt.title('Analytic FFT plot')
        plt.show()

    return sigFFTPos, freqAxisPos


# if __name__ == "__main__":
    
#     df_dados=pd.read_excel('Dados.xlsx')

#     sig = (df_dados['ay'].to_list())
#     FrequenciaDados = 2 # 100 hz

#     dt = 1/FrequenciaDados # 1/dt


#     # Result in frequencies
#     fftPlot(sig, dt=dt)
#     # Result in samples (if the frequencies axis is unknown)
#     fftPlot(sig)
