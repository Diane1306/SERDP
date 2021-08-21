import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def get_wt(wpre, wffp, wpos, tpre, tffp, tpos, method=True):
    wwpre = wpre - wpre.mean()
    if method:
        wwffp = wffp - wpre.mean()
    else:
        wwffp = wffp - wffp.mean()
    wwpos = wpos - wpos.mean()

    ttpre = tpre - tpre.mean()
    if method:
        ttffp = tffp - tpre.mean()
    else:
        ttffp = tffp - tffp.mean()
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
    return Co[:int(n / 2)], Co.sum()


def get_phase(a, b, n):
    afft = np.fft.fft(a, n=n) / len(a)
    bfft = np.fft.fft(b, n=n) / len(b)
    Co = afft.real * bfft.real + afft.imag * bfft.imag
    Qu = afft.imag * bfft.real - afft.real * bfft.imag
    phase = np.arctan2(Qu, Co) * 180. / np.pi + 180.
    # print(phase.max(), phase.min())
    return phase[:int(n / 2)]

def get_quadrant(phi):
    phi = np.array(phi)
    flag1 = np.where(np.logical_or(np.logical_and(phi >= 0, phi < 45), np.logical_and(phi >= 315, phi < 360)), True, False)
    flag2 = np.where(np.logical_and(phi >= 45, phi < 135), True, False)
    flag3 = np.where(np.logical_and(phi >= 135, phi < 225), True, False)
    flag4 = np.where(np.logical_and(phi >= 225, phi < 315), True, False)
    return flag1, flag2, flag3, flag4


def get_smoonthed_cospectra_in_phase(phi, psd):
    flag1, flag2, flag3, flag4 = get_quadrant(phi)
    flag = [flag1, flag2, flag3, flag4]
    # print([f.sum() for f in flag])
    psdr = []
    for fi in range(4):
        psdsmoothed = np.zeros(smoothlen)
        for si in range(smoothlen):
            psdsmoothed[si] = np.nanmean(np.where(flag[fi], psd, np.nan)[smoonthloc[si]:smoonthloc[si + 1]])
        psdr.append(psdsmoothed)
    return psdr


def get_smoothed_slope(psd):
    psdsmoothed = np.zeros(smoothlen)
    for si in range(smoothlen):
        psdsmoothed[si] = psd[smoonthloc[si]:smoonthloc[si + 1]].mean()
    flag = np.where(psdsmoothed > 0, 1, 0)
    psdsmoothed = np.where(flag, psdsmoothed, np.nan)
    slope = 0
    if flag.sum() > 2:
        flagtemp = ~np.isnan(psdsmoothed[28:-4])
        if flagtemp.sum() > 2:
            r = stats.linregress(np.log(FreqSmooth[28:-4][flagtemp]), np.log(psdsmoothed[28:-4][flagtemp]))
            slope = r.slope
    else:
        slope = np.nan
    return psdsmoothed, slope


def get_plot(vn, fn, psd_pre, psd_ffp, psd_pos):
    def plot(ax, yr, f, psd, color):
        flag = ~np.isnan(np.where(psd > 0, psd, np.nan))
        ax.loglog(np.where(psd > 0, f, np.nan)[flag], np.where(psd > 0, psd, np.nan)[flag], color, lw=1)
        ax.set_ylim(10 ** (-8), yr)
        ax.set_xlim(.001, 1.)

    yr = 10
    plt.subplots(1, 3, figsize=(9, 12))
    for hi in range(3):
        ax = plt.subplot(3, 1, hi + 1)
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

        if hi == 2:
            plt.text(.05, .95, '3m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 1:
            plt.text(.05, .95, '10m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 0:
            plt.text(.05, .95, '19m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
    plt.suptitle(f"Cospectrum of {vn}", fontsize='x-large', fontweight='bold', x=.5, y=.98)
    plt.subplots_adjust(top=.95, bottom=.05, right=.95, left=.1,
                        hspace=.1, wspace=0)
    plt.savefig(f'./plot/paper/spectral/EastTower_CoSpectrum{fn}_Averaged_loglog_smoonthed.png', bbox_inches='tight')
    plt.close()


def get_normalized_cospectraplot(vn, fn, phi, psdtot, method):
    def plot(ax, f, psd, color, ls, yl, yr):
        flag = ~np.isnan(psd)
        ax.semilogx(f[flag], psd[flag], color, lw=1, ls=ls)
        ax.set_ylim(yl, yr)
        ax.set_xlim(.001, 1.)
        ax.hlines(0, .001, 1., color='cyan', lw=.8)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)

    color = ['lime', 'r', 'k']
    linestyles = ['-.', ':', '-', '--']
    yl = [-.02, -.01, -.1]
    yr = [.06, .1, .1]
    plt.subplots(3, 3, figsize=(16, 12))
    for hi in range(3):
        for pi in range(3):
            ax = plt.subplot(3, 3, hi * 3 + pi + 1)
            psdsmoothed = get_smoonthed_cospectra_in_phase(phi[hi][pi], psdtot[hi][pi])
            for fi in range(4):
                plot(ax, FreqSmooth, psdsmoothed[fi], color[pi], linestyles[fi], yl[hi], yr[hi])
            if hi == 0:
                plt.xlabel('Frequency ($s^{-1}$)', fontsize=15)
            if pi == 0:
                plt.ylabel("$Co_{w't'}(f)  /  |\overline{w^{'}t^{'}}|$", fontsize=15)


        if hi == 2:
            plt.text(.05, .95, '3m', fontsize=15, fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 1:
            plt.text(.05, .95, '10m', fontsize=15, fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 0:
            plt.text(.05, .95, '19m', fontsize=15, fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        if (hi == 0) and (pi == 2):
            plt.plot([], [], color='k', ls=linestyles[0], label='$0^\circ$')
            plt.plot([], [], color='k', ls=linestyles[1], label='$90^\circ$')
            plt.plot([], [], color='k', ls=linestyles[2], label='$180^\circ$')
            plt.plot([], [], color='k', ls=linestyles[3], label='$270^\circ$')
            plt.legend(loc='upper right', frameon=False, fontsize=15)

        if (hi == 1) and (pi == 2):
            plt.plot([], [], color=color[0], label='Pre-FFP')
            plt.plot([], [], color=color[1], label='FFP')
            plt.plot([], [], color=color[2], label='Pos-FFP')
            plt.legend(loc='upper right', frameon=False, fontsize=15)
    plt.suptitle(f"Cospectrum of {vn}", fontsize=16, fontweight='bold', x=.5, y=.98)
    plt.subplots_adjust(top=.95, bottom=.05, right=.95, left=.1,
                        hspace=.1, wspace=0)
    if method:
        plt.savefig(f'./plot/paper/spectral/EastTower_NormCoSpectrum{fn}_logx_smoonthed.png',
                    bbox_inches='tight')
    else:
        plt.savefig(f'./plot/paper/spectral/PeriodMean_EastTower_NormCoSpectrum{fn}_logx_smoonthed.png',
                    bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    df = pd.read_table('/home/msuclass9/study/Diane/SEDRP/data/FreqForSmooth.txt', delim_whitespace=True,
                       names=('period', 'ps', 'freq', 'cumps', 'sump', 'cumsump', 'index'))
    FreqSmooth = df['freq']
    smoonthloc = [0]
    smoonthloc.extend(df['cumsump'])
    smoothlen = len(FreqSmooth)

    vn = ['Heat Flux total', 'Sweep', 'Ejection', 'Inward Interaction', 'Outward Interaction']
    fn = ['heat_tot', 'Sweep', 'Ejection', 'Inward', 'Outward']

    dfpre = pd.read_csv('/home/msuclass9/study/Diane/SEDRP/data/East_Pre-FFP.csv')
    dfffp = pd.read_csv('/home/msuclass9/study/Diane/SEDRP/data/East_FFP.csv')
    dfpos = pd.read_csv('/home/msuclass9/study/Diane/SEDRP/data/East_Post-FFP.csv')
    method = True
    wwpre3, wwffp3, wwpos3, ttpre3, ttffp3, ttpos3 = get_wt(dfpre['W_3m'], dfffp['W_3m'], dfpos['W_3m'],
                                                            dfpre['T_3m'], dfffp['T_3m'], dfpos['T_3m'], method=method)
    wwpre10, wwffp10, wwpos10, ttpre10, ttffp10, ttpos10 = get_wt(dfpre['W_10m'], dfffp['W_10m'], dfpos['W_10m'],
                                                                  dfpre['T_10m'], dfffp['T_10m'], dfpos['T_10m'],
                                                                  method=method)
    wwpre20, wwffp20, wwpos20, ttpre20, ttffp20, ttpos20 = get_wt(dfpre['W_20m'], dfffp['W_20m'], dfpos['W_20m'],
                                                                  dfpre['T_20m'], dfffp['T_20m'], dfpos['T_20m'],
                                                                  method=method)
    ww = [wwpre20, wwpre10, wwpre3, wwffp20, wwffp10, wwffp3, wwpos20, wwpos10, wwpos3]
    tt = [ttpre20, ttpre10, ttpre3, ttffp20, ttffp10, ttffp3, ttpos20, ttpos10, ttpos3]
    # covar = [(ww[i] * tt[i]).mean() for i in range(9)]
    n = len(wwpre3)
    for ei in range(5):
        Co = []
        phi = []
        for i in range(9):
            wwe, tte = get_event(fn[ei], ww[i], tt[i])
            Cotemp, Covartemp = get_cospectra(wwe, tte, n)
            # Co.append(Cotemp)
            Co.append(Cotemp / Covartemp)
            phi.append(get_phase(wwe, tte, n))
        get_normalized_cospectraplot(vn[ei], fn[ei], [phi[0:3], phi[3:6], phi[6:9]], [Co[0:3], Co[3:6], Co[6:9]],
                                     method)
        # get_plot(vn[ei], fn[ei], Co[0:3], Co[3:6], Co[6:9])
