import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def get_wt(wpre, wffp, wpos, tpre, tffp, tpos):
    wwpre = wpre - wpre.mean()
    wwffp = wffp - wpre.mean()
    wwpos = wpos - wpos.mean()

    ttpre = tpre - tpre.mean()
    ttffp = tffp - tpre.mean()
    ttpos = tpos - tpos.mean()
    return wwpre, wwffp, wwpos, ttpre, ttffp, ttpos


def get_event(en, ww, tt):
    if en == 'heat_tot':
        wwr = ww
        ttr = tt
    elif en == 'Sweep':
        wwr = np.where(np.logical_and(ww<0, tt<0), ww, 0)
        ttr = np.where(np.logical_and(ww<0, tt<0), tt, 0)
    elif en == 'Ejection':
        wwr = np.where(np.logical_and(ww>0, tt>0), ww, 0)
        ttr = np.where(np.logical_and(ww>0, tt>0), tt, 0)
    elif en == 'Inward':
        wwr = np.where(np.logical_and(ww<0, tt>0), ww, 0)
        ttr = np.where(np.logical_and(ww<0, tt>0), tt, 0)
    elif en == 'Outward':
        wwr = np.where(np.logical_and(ww>0, tt<0), ww, 0)
        ttr = np.where(np.logical_and(ww>0, tt<0), tt, 0)
    return wwr, ttr


def get_cospectra(a, b):
    afft = np.fft.fft(a) / len(a)
    bfft = np.fft.fft(b) / len(b)
    Co = afft.real * bfft.real + afft.imag * bfft.imag
    Qu = afft.imag * bfft.real - afft.real * bfft.imag
    return Co, Qu




