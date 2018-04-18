#!/usr/bin/env python
# -*- coding: utf-8 -*-


_author__ = "Patrice Guyot"
__email__ = "patrice.guyot@irit.fr"
__version__ = "1.1.0"

"""
***********************************************************************
Name            : audio_processing.py
Description     : The module implements a function to read audio file.
Author          : Patrice Guyot.
***********************************************************************
"""



from scipy.io.wavfile import read as wavread
import logging
from os import remove, sep
import subprocess



def audio_read(audioFilePath, formatsox=False):
    """

    Read an audio file (from scipy.io.wavfile)

    A conversation of the aufio file can be processed using SOX (by default set at False).

    :param audioFilePath: audio file name (with eventually the full path)
    :type audioFilePath: str
    :param formatsox: if set at True, us SOX to convert audio file in "wav file, 1 channel, 16 kHz, 16bits". Default: False
    :type audioFilePath: str

    :returns:
        * sr: sampling rate of the signal
        * sig: list of values of the signal
    :rtype: tuple

    """
    logging.info('Reading of the audio file : ' + audioFilePath)
    if formatsox:
        tmpFile = "tmp.wav"
        logging.info('\t- Conversion en wav 16k with SOX.')
        cmd = "sox " + audioFilePath + " -c 1 -r 16k -b 16 -G " + tmpFile + " rate -m"
        subprocess.check_call(cmd)
        [sr, sig] = wavread(tmpFile)
        remove(tmpFile)
    else:
        [sr, sig] = wavread(audioFilePath)

    return sr, sig

