import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_data(contri, fracti):
    EventMean = np.zeros((3, 5, 3, 4))
    Frac = np.zeros((3, 5, 3, 4))
    for hi in range(3):
        for ti in range(5):
            for pi in range(3):
                for ei in range(4):
                    if ei < 2:
                        EventMean[hi, ti, pi, ei] = contri[ti, ei, hi * 3 + pi]
                        Frac[hi, ti, pi, ei] = fracti[ti, ei, hi * 3 + pi]
                    elif ei == 2:
                        EventMean[hi, ti, pi, ei] = contri[ti, ei + 1, hi * 3 + pi]
                        Frac[hi, ti, pi, ei] = fracti[ti, ei + 1, hi * 3 + pi]
                    elif ei == 3:
                        EventMean[hi, ti, pi, ei] = contri[ti, ei - 1, hi * 3 + pi]
                        Frac[hi, ti, pi, ei] = fracti[ti, ei - 1, hi * 3 + pi]

    EventMean = np.nanmean(EventMean, axis=1)
    HeatMean = EventMean.sum(axis=2)
    Frac = np.nanmean(Frac, axis=1)
    return EventMean, HeatMean, Frac

def get_plot(EventMeanc, HeatMeanc, Fracc, EventMeanf, HeatMeanf, Fracf):
    # height(3,10,20), period(Pre, FFP, Post), event(Sweep, Ejection, Inward, Outward)
    fig, axs = plt.subplots(2, 3, figsize=(9, 5))
    plt.subplots_adjust(wspace=0.05, hspace=0.35, bottom=0.0)
    Y = [10, 19]
    period = ['Pre-FFP', 'FFP', 'Post-FFP']
    for pi in range(3):
        ax = plt.subplot(2, 3, pi + 1)
        plt.barh(Y, EventMeanc[1:, pi, 0] + EventMeanc[1:, pi, 1], align='center', height=2, color='pink')
        plt.barh(Y, EventMeanc[1:, pi, 0], align='center', height=2, color='red')
        plt.barh(Y, EventMeanc[1:, pi, 2] + EventMeanc[1:, pi, 3], align='center', height=2, color='b')
        plt.barh(Y, EventMeanc[1:, pi, 2], align='center', height=2, color='cyan')

        plt.vlines(HeatMeanc[1, pi], 8.5, 11.5, lw=1, colors='yellow')
        plt.text(HeatMeanc[1, pi] - .018, 11.7, f'{HeatMeanc[1, pi]:.3f}', fontsize=10)
        plt.vlines(HeatMeanc[2, pi], 17.5, 20.5, lw=1, colors='yellow')
        plt.text(HeatMeanc[2, pi] - .018, 20.7, f'{HeatMeanc[2, pi]:.3f}', fontsize=10)
        plt.vlines(0, 0, 22, lw=1, ls='dashed', color='grey')
        plt.xlim(-.4, .7)
        plt.xticks([i * .1 for i in range(-4, 7)], ['-0.4', '', '-0.2', '', '0.0', '', '0.2', '', '0.4', '', '0.6', ''],
                   fontsize=11)
        plt.ylim(0, 22)
        if pi == 0:
            plt.ylabel("Height (m)", fontsize=12)
            plt.yticks([i for i in range(0, 22, 3)], fontsize=11)
            plt.text(-.2, 1, '(a)', fontsize=15, fontweight='bold', c='k', ha='center', va='center',
                     transform=ax.transAxes)
        else:
            plt.yticks([i for i in range(0, 22, 3)], [])


        axx = ax.twiny()
        axx.plot(Fracc[1:, pi, 0], Y, 'r-o', ms=3)
        axx.scatter(Fracc[1:, pi, 0], Y, color='m')
        axx.plot(Fracc[1:, pi, 1], Y, c='pink', ls='-', marker='o', ms=3)
        axx.scatter(Fracc[1:, pi, 1], Y, color='m')
        axx.plot(Fracc[1:, pi, 2], Y, c='cyan', ls='-', marker='o', ms=3)
        axx.scatter(Fracc[1:, pi, 2], Y, color='m')
        axx.plot(Fracc[1:, pi, 3], Y, 'b-o', ms=3)
        axx.scatter(Fracc[1:, pi, 3], Y, color='m')
        axx.set_xlim(0, .7)
        axx.set_xticks(np.arange(0, 0.7, 0.1))
        axx.set_xticklabels(['0.0', '', '0.2', '', '0.4', '', '0.6', ''], fontsize=11)
        axx.set_xlabel(f'Fraction of Events ({period[pi]})', fontsize=12)
        axx.set_yticks([3, 10, 19])
        plt.legend(['Sweep', 'Ejection', 'Inward', 'Outward'], loc='lower right', frameon=False, fontsize=9)

    Y = [3, 10, 19]
    for pi in range(3):
        ax = plt.subplot(2, 3, 3+ pi + 1)
        plt.barh(Y, EventMeanf[:, pi, 0] + EventMeanf[:, pi, 1], align='center', height=2, color='pink')
        plt.barh(Y, EventMeanf[:, pi, 0], align='center', height=2, color='red')
        plt.barh(Y, EventMeanf[:, pi, 2] + EventMeanf[:, pi, 3], align='center', height=2, color='b')
        plt.barh(Y, EventMeanf[:, pi, 2], align='center', height=2, color='cyan')

        plt.vlines(HeatMeanf[0, pi], 1.5, 4.5, lw=1, colors='yellow')
        plt.text(HeatMeanf[0, pi] - .018, 4.7, f'{HeatMeanf[0, pi]:.3f}', fontsize=10)
        plt.vlines(HeatMeanf[1, pi], 8.5, 11.5, lw=1, colors='yellow')
        plt.text(HeatMeanf[1, pi] - .018, 11.7, f'{HeatMeanf[1, pi]:.3f}', fontsize=10)
        plt.vlines(HeatMeanf[2, pi], 17.5, 20.5, lw=1, colors='yellow')
        plt.text(HeatMeanf[2, pi] - .018, 20.7, f'{HeatMeanf[2, pi]:.3f}', fontsize=10)

        plt.vlines(0, 0, 22, lw=1, ls='dashed', color='grey')
        if pi-1:
            plt.xlim(-.4, .7)
            plt.xticks([i * .1 for i in range(-4, 7)], ['-0.4', '', '-0.2', '', '0.0', '', '0.2', '', '0.4', '', '0.6', ''],
                       fontsize=11)
        else:
            plt.xlim(-1, 5)
            plt.xticks([i for i in range(-1, 5)], fontsize=11)
        plt.xlabel("$w't' (\u2103ms^{-1})$ " +f"({period[pi]})", fontsize=12)
        plt.ylim(0, 22)
        if pi == 0:
            plt.ylabel("Height (m)", fontsize=12)
            plt.yticks([i for i in range(0, 22, 5)], fontsize=11)
            plt.text(-.2, 1, '(b)', fontsize=15, fontweight='bold', c='k', ha='center', va='center',
                     transform=ax.transAxes)
        else:
            plt.yticks([i for i in range(0, 22, 5)], [])


        axx = ax.twiny()
        axx.plot(Fracf[:, pi, 0], Y, 'r-o', ms=3)
        axx.scatter(Fracf[:, pi, 0], Y, color='m')
        axx.plot(Fracf[:, pi, 1], Y, c='pink', ls='-', marker='o', ms=3)
        axx.scatter(Fracf[:, pi, 1], Y, color='m')
        axx.plot(Fracf[:, pi, 2], Y, c='cyan', ls='-', marker='o', ms=3)
        axx.scatter(Fracf[:, pi, 2], Y, color='m')
        axx.plot(Fracf[:, pi, 3], Y, 'b-o', ms=3)
        axx.scatter(Fracf[:, pi, 3], Y, color='m')
        axx.set_xlim(0, .7)
        axx.set_xticks(np.arange(0, 0.7, 0.1))
        axx.set_xticklabels(['0.0', '', '0.2', '', '0.4', '', '0.6', ''], fontsize=11)
        axx.set_yticks([3, 10, 19])

    plt.savefig('../../2020_SERDP_Diane/plot/ContributionFreq_AverControlFluxTower.png', bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    contri = np.load('../../2020_SERDP_Diane/data/ContributionValue_Control.npy')  # (5, 4, 9)
    fracti = np.load('../../2020_SERDP_Diane/data/FractionValue_Control.npy')
    EventMeanc, HeatMeanc, Fracc = get_data(contri, fracti)

    contri = np.load('../../2020_SERDP_Diane/data/ContributionValue_Flux.npy')  # (5, 4, 9)
    fracti = np.load('../../2020_SERDP_Diane/data/FractionValue_Flux.npy')
    EventMeanf, HeatMeanf, Fracf = get_data(contri, fracti)
    get_plot(EventMeanc, HeatMeanc, Fracc, EventMeanf, HeatMeanf, Fracf)


