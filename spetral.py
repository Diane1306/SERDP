import os

import pandas as pd
from scipy import stats
import numpy as np
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt

def get_data(dr, height, tower):
    DIR = f'{dr}{height}'
    fna = [name for name in os.listdir(DIR)]
    fna.sort()

    ww = []
    tt = []
    Sweep = []
    Ejection = []
    Outward = []
    Inward = []
    CumSec = []

    for fi in range(len(tower)):
        xls = pd.ExcelFile(f'{dr}{height}/{fna[fi]}')
        df_pre = pd.read_excel(xls, 'Pre-FFP')
        df_dur = pd.read_excel(xls, 'FFP')
        df_post = pd.read_excel(xls, 'Post-FFP')

        df = [df_pre, df_dur, df_post]
        for di in range(3):
            ww.append(df[di]["w' (m/s)"].dropna())
            tt.append(df[di]["t' (C)"].dropna())
            Sweep.append(df[di]["w't' (Sweep)"].dropna())
            Ejection.append(df[di]["w't' (Ejection)"].dropna())
            Outward.append(df[di]["w't' (Out Int)"].dropna())
            Inward.append(df[di]["w't' (In Int)"].dropna())
            CumSec.append(df[di].iloc[:, 0].dropna())
    return ww, tt, Sweep, Ejection, Outward, Inward, CumSec


def calc_slope1(psd):
    # slope for each tower
    psdsmoonthed = np.zeros(smoothlen)
    for si in range(smoothlen):
        psdsmoonthed[si] = psd[smoonthloc[si]:smoonthloc[si + 1]].mean()
    r = stats.linregress(np.log(FreqSmooth[28:-4]), np.log(psdsmoonthed[28:-4]))
    return r.slope

def calc_slope(psd):
    # slope for averaged psd
    r = stats.linregress(np.log(FreqSmooth[28:-4]), np.log(psd[28:-4]))
    return r.slope


def calc_LombScargleSpectrum(cumsec, flux):
    psd_pre = 0
    psd_ffp = 0
    psd_pos = 0
    slope = np.zeros(3)
    tn = int(len(flux) / 3)
    for ti in range(tn):
        psdtemp = LombScargle(cumsec[ti * 3 + 0], flux[0 + 3 * ti], normalization='psd').power(freq)
        psd_pre += psdtemp
        slope[0] += calc_slope1(psdtemp)
        psdtemp = LombScargle(cumsec[ti * 3 + 1], flux[1 + 3 * ti], normalization='psd').power(freq)
        psd_ffp += psdtemp
        slope[1] += calc_slope1(psdtemp)
        psdtemp = LombScargle(cumsec[ti * 3 + 2], flux[2 + 3 * ti], normalization='psd').power(freq)
        psd_pos += psdtemp
        slope[2] += calc_slope1(psdtemp)
    psd_pre = psd_pre / tn
    psd_ffp = psd_ffp / tn
    psd_pos = psd_pos / tn
    slope[0] = slope[0] / tn
    slope[1] = slope[1] / tn
    slope[2] = slope[2] / tn
    return psd_pre, psd_ffp, psd_pos, slope[0], slope[1], slope[2]


def get_plot(vn, fn, psd, variance):
    psd = [psd[i] / variance[i] for i in range(9)]
    ylim = .1
    def plot(ax, yr, f, psd, color):
        ax.semilogx(f, psd, color, lw=1)
        ax.set_ylim(0, yr)
        ax.set_xlim(.001, 0.5)

    yr = [ylim for i in range(3)]
    plt.subplots(1, 3, figsize=(9, 12))
    for hi in range(3):
        ax = plt.subplot(3, 1, 3 - hi)
        plot(ax, yr[hi], freq, psd[hi * 3], 'lime')
        plot(ax, yr[hi], freq, psd[1 + hi * 3], 'r')
        plot(ax, yr[hi], freq, psd[2 + hi * 3], 'k')
        if hi == 0:
            plt.xlabel('Frequency ($s^{-1}$)', fontsize='large')
        plt.ylabel('norm PSD', fontsize='large')
        plt.legend([f'Pre-FFP ({variance[hi * 3] / len(freq):.2f} $\u2103^{2}m^{2}s^{{-2}})$',
                    f'FFP ({variance[hi * 3 + 1] / len(freq):.2f} $\u2103^{2}m^{2}s^{{-2}})$',
                    f'Post-FFP ({variance[hi * 3 + 2] / len(freq):.2f} $\u2103^{2}m^{2}s^{{-2}})$'], loc='upper right',
                   frameon=False, fontsize='large')

        if hi == 0:
            plt.text(.05, .95, '3m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 1:
            plt.text(.05, .95, '10m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
        elif hi == 2:
            plt.text(.05, .95, '19m', fontsize='large', fontweight='bold', ha='center', va='center',
                     transform=ax.transAxes)
    plt.suptitle(f"Normalized Power Spectrum Density of {vn}", fontsize='x-large', fontweight='bold', x=.5, y=.98)
    plt.subplots_adjust(top=.95, bottom=.05, right=.95, left=.1,
                        hspace=.1, wspace=0)
    plt.savefig(f'./plot/paper/spectral/PSD{fn}_Averaged_logx.png', bbox_inches='tight')
    plt.close()
    return

if __name__ == "__main__":
    data_dir = '/home/msuclass9/study/Diane/SEDRP/data/PeriodMean/'
    ww20, tt20, Sweep20, Ejection20, Outward20, Inward20, CumSec20 = get_data(data_dir, '20m',
                                                                              ['East', 'Flux', 'North', 'South',
                                                                               'West'])
    ww10, tt10, Sweep10, Ejection10, Outward10, Inward10, CumSec10 = get_data(data_dir, '10m',
                                                                              ['East', 'Flux', 'North', 'West'])
    ww3, tt3, Sweep3, Ejection3, Outward3, Inward3, CumSec3 = get_data(data_dir, '3m',
                                                                       ['East', 'Flux', 'North', 'South',
                                                                        'West'])
    # heat flux total
    heat_tot3 = [Sweep3[i].values + Ejection3[i].values + Inward3[i].values + Outward3[i].values for i in range(15)]
    heat_tot10 = [Sweep10[i].values + Ejection10[i].values + Inward10[i].values + Outward10[i].values for i in
                  range(12)]
    heat_tot20 = [Sweep20[i].values + Ejection20[i].values + Inward20[i].values + Outward20[i].values for i in
                  range(15)]
    flux3 = [heat_tot3, Sweep3, Ejection3, Inward3, Outward3]
    flux10 = [heat_tot10, Sweep10, Ejection10, Inward10, Outward10]
    flux20 = [heat_tot20, Sweep20, Ejection20, Inward20, Outward20]

    # get freq for smooth
    df = pd.read_table('./data/FreqForSmooth.txt', delim_whitespace=True,
                       names=('period', 'ps', 'freq', 'cumps', 'sump', 'cumsump', 'index'))
    FreqSmooth = df['freq']
    smoonthloc = [0]
    smoonthloc.extend(df['cumsump'])
    smoothlen = len(FreqSmooth)

    freq = np.linspace(0, 5, 9001)[1:]
    vn = ['Heat Flux total', 'Sweep', 'Ejection', 'Inward Interaction', 'Outward Interaction']
    fn = ['heat_tot', 'Sweep', 'Ejection', 'Inward', 'Outward']
    for ei in range(5):
        slope = np.zeros(9)
        psd3_pre, psd3_ffp, psd3_pos, slope[0], slope[1], slope[2] = calc_LombScargleSpectrum(CumSec3, flux3[ei])
        psd10_pre, psd10_ffp, psd10_pos, slope[3], slope[4], slope[5] = calc_LombScargleSpectrum(CumSec10, flux10[ei])
        psd20_pre, psd20_ffp, psd20_pos, slope[6], slope[7], slope[8] = calc_LombScargleSpectrum(CumSec20, flux20[ei])

        variance = [np.nansum(psd3_pre), np.nansum(psd3_ffp), np.nansum(psd3_pos), np.nansum(psd10_pre),
                    np.nansum(psd10_ffp),
                    np.nansum(psd10_pos), np.nansum(psd20_pre), np.nansum(psd20_ffp), np.nansum(psd20_pos)]

        psd = [psd3_pre, psd3_ffp, psd3_pos, psd10_pre, psd10_ffp, psd10_pos, psd20_pre, psd20_ffp, psd20_pos]
        get_plot(vn[ei], fn[ei], psd, variance)