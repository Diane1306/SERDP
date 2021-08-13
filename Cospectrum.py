import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


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
        wwr = np.where(np.logical_and(ww < 0, tt < 0), ww, 0)
        ttr = np.where(np.logical_and(ww < 0, tt < 0), tt, 0)
    elif en == 'Ejection':
        wwr = np.where(np.logical_and(ww > 0, tt > 0), ww, 0)
        ttr = np.where(np.logical_and(ww > 0, tt > 0), tt, 0)
    elif en == 'Inward':
        wwr = np.where(np.logical_and(ww < 0, tt > 0), ww, 0)
        ttr = np.where(np.logical_and(ww < 0, tt > 0), tt, 0)
    elif en == 'Outward':
        wwr = np.where(np.logical_and(ww > 0, tt < 0), ww, 0)
        ttr = np.where(np.logical_and(ww > 0, tt < 0), tt, 0)
    return wwr, ttr


def get_cospectra(a, b, n):
    afft = np.fft.fft(a, n=n) / len(a)
    bfft = np.fft.fft(b, n=n) / len(b)
    Co = afft.real * bfft.real + afft.imag * bfft.imag
    # Qu = afft.imag * bfft.real - afft.real * bfft.imag
    return Co[:int(n / 2)]


def get_smoothed_slope(psd):
    df = pd.read_table('./data/FreqForSmooth.txt', delim_whitespace=True,
                       names=('period', 'ps', 'freq', 'cumps', 'sump', 'cumsump', 'index'))
    freqsmooth = df['freq']
    smoonthloc = [0]
    smoonthloc.extend(df['cumsump'])
    smoothlen = len(freqsmooth)
    psdsmoothed = np.zeros(smoothlen)
    for si in range(smoothlen):
        psdsmoothed[si] = psd[smoonthloc[si]:smoonthloc[si + 1]].mean()
    flag = np.where(psdsmoothed > 0, 1, 0)
    psdsmoothed = np.where(flag, psdsmoothed, np.nan)
    slope = 0
    if flag.sum() > 2:
        flagtemp = ~np.isnan(psdsmoothed[28:-4])
        if flagtemp.sum() > 2:
            r = stats.linregress(np.log(freqsmooth[28:-4][flagtemp]), np.log(psdsmoothed[28:-4][flagtemp]))
            slope = r.slope
    else:
        slope = np.nan
    return psdsmoothed, slope


def get_plot(vn, fn, psd_pre, psd_ffp, psd_pos):
    def plot(ax, yr, f, psd, color):
        ax.loglog(np.where(psd > 0, f, np.nan), np.where(psd > 0, psd, np.nan), color, lw=1)
        ax.set_ylim(10 ** (-8), yr)
        ax.set_xlim(.001, 1.)

    yr = 10
    plt.subplots(1, 3, figsize=(9, 12))
    for hi in range(3):
        ax = plt.subplot(3, 1, 3 - hi)
        psdpre, slopepre = get_smoothed_slope(psd_pre[hi])
        psdffp, slopeffp = get_smoothed_slope(psd_ffp[hi])
        psdpos, slopepos = get_smoothed_slope(psd_pos[hi])
        plot(ax, yr, FreqSmooth, psdpre, 'lime')
        plot(ax, yr, FreqSmooth, psdffp, 'r')
        plot(ax, yr, FreqSmooth, psdpos, 'k')
        ax.loglog(FreqSmooth[28:], FreqSmooth[28:] ** (-5 / 3) * psdffp[28] / 50, 'm--')
        ax.loglog(FreqSmooth[28:], FreqSmooth[28:] ** (-7 / 3) * psdffp[28] / 1000, 'c--')
        if hi == 0:
            plt.xlabel('Frequency ($s^{-1}$)', fontsize='large')
        plt.ylabel('Cospectrum ($\u2103^{2}m^{2}s^{{-2}}$)', fontsize='large')
        plt.legend([f'Pre-FFP ({slopepre:.2f})', f'FFP ({slopeffp:.2f})', f'Post-FFP ({slopepos:.2f})',
                    '$k^{-5/3}$' + f' ({-5 / 3:.2f})', '$k^{-7/3}$' + f' ({-7 / 3:.2f})'],
                   loc='upper right', frameon=False, fontsize='large')

        if hi == 0:
            plt.text(.05, .95, '3m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 1:
            plt.text(.05, .95, '10m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 2:
            plt.text(.05, .95, '19m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
    plt.suptitle(f"Cospectrum of {vn}", fontsize='x-large', fontweight='bold', x=.5, y=.98)
    plt.subplots_adjust(top=.95, bottom=.05, right=.95, left=.1,
                        hspace=.1, wspace=0)
    plt.savefig(f'./plot/paper/spectral/EastTower_CoSpectrum{fn}_Averaged_loglog_smoonthed.png', bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    df = pd.read_table('./data/FreqForSmooth.txt', delim_whitespace=True,
                       names=('period', 'ps', 'freq', 'cumps', 'sump', 'cumsump', 'index'))
    FreqSmooth = df['freq']
    smoothlen = len(FreqSmooth)
    vn = ['Heat Flux total', 'Sweep', 'Ejection', 'Inward Interaction', 'Outward Interaction']
    fn = ['heat_tot', 'Sweep', 'Ejection', 'Inward', 'Outward']

    dfpre = pd.read_csv('./data/East_Pre-FFP.csv')
    dfffp = pd.read_csv('./data/East_FFP.csv')
    dfpos = pd.read_csv('./data/East_Post-FFP.csv')
    wwpre3, wwffp3, wwpos3, ttpre3, ttffp3, ttpos3 = get_wt(dfpre['W_3m'], dfffp['W_3m'], dfpos['W_3m'],
                                                            dfpre['T_3m'], dfffp['T_3m'], dfpos['T_3m'])
    wwpre10, wwffp10, wwpos10, ttpre10, ttffp10, ttpos10 = get_wt(dfpre['W_10m'], dfffp['W_10m'], dfpos['W_10m'],
                                                                  dfpre['T_10m'], dfffp['T_10m'], dfpos['T_10m'])
    wwpre20, wwffp20, wwpos20, ttpre20, ttffp20, ttpos20 = get_wt(dfpre['W_20m'], dfffp['W_20m'], dfpos['W_20m'],
                                                                  dfpre['T_20m'], dfffp['T_20m'], dfpos['T_20m'])
    ww = [wwpre3, wwpre10, wwpre20, wwffp3, wwffp10, wwffp20, wwpos3, wwpos10, wwpos20]
    tt = [ttpre3, ttpre10, ttpre20, ttffp3, ttffp10, ttffp20, ttpos3, ttpos10, ttpos20]
    n = len(wwpre3)
    for ei in range(5):
        Co = []
        for i in range(9):
            wwe, tte = get_event(fn[ei], ww[i], tt[i])
            Co.append(get_cospectra(wwe, tte, n))
        get_plot(vn[ei], fn[ei], Co[0:3], Co[3:6], Co[6:9])
