#!/usr/bin/env python
# -*- coding: utf-8 -*-



_author__ = "Patrice Guyot"
__email__ = "patrice.guyot@irit.fr"
__version__ = "1.1.0"


"""
***********************************************************************
Name            : yin.py
Description     : Fundamental frequency estimation. Based on the YIN alorgorithm [1]: De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.
Author          : Patrice Guyot. Previous works on the implementation of the YIN algorithm have been made at IRIT thanks to Robin Larvor, Maxime Le Coz and Lionel Koenig.
***********************************************************************

"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from os import sep
from audio_processing import *
import time






def differenceFunction_original(x, N, tau_max):
    """
    Compute difference function of data x. This corresponds to equation (6) in [1]

    Original algorithm.

    :param x: audio data
    :param N: length of data
    :param tau_max: integration window size
    :return: difference function
    :rtype: list
    """
    df = [0] * tau_max
    for tau in range(1, tau_max):
         for j in range(0, N - tau_max):
             tmp = long(x[j] - x[j + tau])
             df[tau] += tmp * tmp
    return df


def differenceFunction_scipy(x, N, tau_max):
    """
    Compute difference function of data x. This corresponds to equation (6) in [1]

    Faster implementation of the difference function.
    The required calculation can be easily evaluated by Autocorrelation function or similarly by convolution.
    Wiener–Khinchin theorem allows computing the autocorrelation with two Fast Fourier transforms (FFT), with time complexity O(n log n).
    This function use an accellerated convolution function fftconvolve from Scipy package.

    :param x: audio data
    :param N: length of data
    :param tau_max: integration window size
    :return: difference function
    :rtype: list
    """
    x = np.array(x, np.float64)
    w = x.size
    x_cumsum = np.concatenate((np.array([0]), (x * x).cumsum()))
    conv = fftconvolve(x, x[::-1])
    tmp = x_cumsum[w:0:-1] + x_cumsum[w] - x_cumsum[:w] - 2 * conv[w - 1:]
    return tmp[:tau_max + 1]




def differenceFunction(x, N, tau_max):
    """
    Compute difference function of data x. This corresponds to equation (6) in [1]

    Fastest implementation. Use the same approach than differenceFunction_scipy.
    This solution is implemented directly with Numpy fft.


    :param x: audio data
    :param N: length of data
    :param tau_max: integration window size
    :return: difference function
    :rtype: list
    """

    x = np.array(x, np.float64)
    w = x.size
    tau_max = min(tau_max, w)
    x_cumsum = np.concatenate((np.array([0.]), (x * x).cumsum()))
    size = w + tau_max
    p2 = (size // 32).bit_length()
    nice_numbers = (16, 18, 20, 24, 25, 27, 30, 32)
    size_pad = min(x * 2 ** p2 for x in nice_numbers if x * 2 ** p2 >= size)
    fc = np.fft.rfft(x, size_pad)
    conv = np.fft.irfft(fc * fc.conjugate())[:tau_max]
    return x_cumsum[w:w - tau_max:-1] + x_cumsum[w] - x_cumsum[:tau_max] - 2 * conv



def cumulativeMeanNormalizedDifferenceFunction(df, N):
    """
    Compute cumulative mean normalized difference function (CMND).

    This corresponds to equation (8) in [1]

    :param df: Difference function
    :param N: length of data
    :return: cumulative mean normalized difference function
    :rtype: list
    """

    cmndf = df[1:] * range(1, N) / np.cumsum(df[1:]).astype(float) #scipy method
    return np.insert(cmndf, 0, 1)



def getPitch(cmdf, tau_min, tau_max, harmo_th=0.1):
    """
    Return fundamental period of a frame based on CMND function.

    :param cmdf: Cumulative Mean Normalized Difference function
    :param tau_min: minimum period for speech
    :param tau_max: maximum period for speech
    :param harmo_th: harmonicity threshold to determine if it is necessary to compute pitch frequency
    :return: fundamental period if there is values under threshold, 0 otherwise
    :rtype: float
    """
    tau = tau_min
    while tau < tau_max:
        if cmdf[tau] < harmo_th:
            while tau + 1 < tau_max and cmdf[tau + 1] < cmdf[tau]:
                tau += 1
            return tau
        tau += 1

    return 0    # if unvoiced



def compute_yin(sig, sr, dataFileName=None, w_len=512, w_step=256, f0_min=100, f0_max=500, harmo_thresh=0.1):
    """

    Compute the Yin Algorithm. Return fundamental frequency and harmonic rate.

    :param sig: Audio signal (list of float)
    :param sr: sampling rate (int)
    :param w_len: size of the analysis window (samples)
    :param w_step: size of the lag between two consecutives windows (samples)
    :param f0_min: Minimum fundamental frequency that can be detected (hertz)
    :param f0_max: Maximum fundamental frequency that can be detected (hertz)
    :param harmo_tresh: Threshold of detection. The yalgorithmù return the first minimum of the CMND fubction below this treshold.

    :returns:

        * pitches: list of fundamental frequencies,
        * harmonic_rates: list of harmonic rate values for each fundamental frequency value (= confidence value)
        * argmins: minimums of the Cumulative Mean Normalized DifferenceFunction
        * times: list of time of each estimation
    :rtype: tuple
    """

    print('\t- Yin: compute yin algorithm')
    tau_min = int(sr / f0_max)
    tau_max = int(sr / f0_min)

    timeScale = range(0, len(sig) - w_len, w_step)  # time values for each analysis window
    times = [t/float(sr) for t in timeScale]
    frames = [sig[t:t + w_len] for t in timeScale]

    pitches = [0.0] * len(timeScale)
    harmonic_rates = [0.0] * len(timeScale)
    argmins = [0.0] * len(timeScale)

    for i, frame in enumerate(frames):

        #Compute YIN
        df = differenceFunction(frame, w_len, tau_max)
        cmdf = cumulativeMeanNormalizedDifferenceFunction(df, tau_max)
        p = getPitch(cmdf, tau_min, tau_max, harmo_thresh)

        #Get results
        if np.argmin(cmdf)>tau_min:
            argmins[i] = float(sr / np.argmin(cmdf))
        if p != 0: # A pitch was found
            pitches[i] = float(sr / p)
            harmonic_rates[i] = cmdf[p]
        else: # No pitch, but we compute a value of the harmonic rate
            harmonic_rates[i] = min(cmdf)


    if dataFileName is not None:
        np.savez(dataFileName, times=times, sr=sr, w_len=w_len, w_step=w_step, f0_min=f0_min, f0_max=f0_max, harmo_thresh=harmo_thresh, pitches=pitches, harmonic_rates=harmonic_rates, argmins=argmins)
        print('\t\t- Data file written in: ' + dataFileName)

    return pitches, harmonic_rates, argmins, times


def main(audioFileName="out.wav", w_len=1024, w_step=256, f0_min=70, f0_max=200, harmo_thresh=0.85, audioDir="./", dataFileName=None, verbose=4):
    """
    Run the computation of the Yin algorithm on a example file.

    Write the results (pitches, harmonic rates, parameters ) in a numpy file.

    :param audioFileName: name of the audio file
    :type audioFileName: str
    :param w_len: length of the window
    :type wLen: int
    :param wStep: length of the "hop" size
    :type wStep: int
    :param f0_min: minimum f0 in Hertz
    :type f0_min: float
    :param f0_max: maximum f0 in Hertz
    :type f0_max: float
    :param harmo_thresh: harmonic threshold
    :type harmo_thresh: float
    :param audioDir: path of the directory containing the audio file
    :type audioDir: str
    :param dataFileName: file name to output results
    :type dataFileName: str
    :param verbose: Outputs on the console : 0-> nothing, 1-> warning, 2 -> info, 3-> debug(all info), 4 -> plot + all info
    :type verbose: int
    """

    if audioDir is not None:
        audioFilePath = audioDir + sep + audioFileName
    else:
        audioFilePath = audioFileName

    sr, sig = audio_read(audioFilePath, formatsox=False)

    start = time.time()
    pitches, harmonic_rates, argmins, times = compute_yin(sig, sr, dataFileName, w_len, w_step, f0_min, f0_max, harmo_thresh)
    end = time.time()
    print("Yin computed in : ", end - start)

    duration = len(sig)/float(sr)

    if verbose >3:
        ax1 = plt.subplot(4, 1, 1)
        ax1.plot([float(x) * duration / len(sig) for x in range(0, len(sig))], sig)
        ax1.set_title('Audio data')
        ax1.set_ylabel('Amplitude')
        ax2 = plt.subplot(4, 1, 2)
        ax2.plot([float(x) * duration / len(pitches) for x in range(0, len(pitches))], pitches)
        ax2.set_title('F0')
        ax2.set_ylabel('Frequency (Hz)')
        ax3 = plt.subplot(4, 1, 3, sharex=ax2)
        ax3.plot([float(x) * duration / len(harmonic_rates) for x in range(0, len(harmonic_rates))], harmonic_rates)
        ax3.plot([float(x) * duration / len(harmonic_rates) for x in range(0, len(harmonic_rates))], [harmo_thresh] * len(harmonic_rates), 'r')
        ax3.set_title('Harmonic rate')
        ax3.set_ylabel('Rate')
        ax4 = plt.subplot(4, 1, 4, sharex=ax2)
        ax4.plot([float(x) * duration / len(argmins) for x in range(0, len(argmins))], argmins)
        ax4.set_title('Index of minimums of CMND')
        ax4.set_ylabel('Frequency (Hz)')
        ax4.set_xlabel('Time (seconds)')
        plt.show()



if __name__ == '__main__':
    main()


