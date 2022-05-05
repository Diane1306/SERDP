import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_table():
    df3 = pd.read_table('~/study/2020_SERDP_Diane/data/TWR_C_2nd_cook_03m_fft1.txt', delim_whitespace=True, index_col=0,
                        names=['Frac', 'U', 'V', 'W', 'T', 'year', 'mon', 'day', 'hour', 'min', 'tenth'])
    df10 = pd.read_table('~/study/2020_SERDP_Diane/data/TWR_C_2nd_cook_10m_fft1.txt', delim_whitespace=True, index_col=0,
                         names=['Frac', 'U', 'V', 'W', 'T', 'year', 'mon', 'day', 'hour', 'min', 'tenth'])
    df20 = pd.read_table('~/study/2020_SERDP_Diane/data/TWR_C_2nd_cook_20m_fft1.txt', delim_whitespace=True, index_col=0,
                         names=['Frac', 'U', 'V', 'W', 'T', 'year', 'mon', 'day', 'hour', 'min', 'tenth'])
    return df3, df10, df20


def get_data(df3, df10, df20):
    # period mean for pre and post and pre mean for ffp
    t = []
    w = []
    wt = []
    tmean = []
    wmean = []
    wtmean = []

    df = [df3, df10, df20]
    for hi in range(3):
        df_pre = df[hi].iloc[cstart[ti]:cffps[ti]]  # pre
        df_ffp = df[hi].iloc[cffps[ti]:cffpe[ti]]  # ffp
        df_pos = df[hi].iloc[cffpe[ti]:cend[ti]]  # post

        t.append(df_pre['T'] - df_pre['T'].mean())
        t.append(df_ffp['T'] - df_pre['T'].mean())
        t.append(df_pos['T'] - df_pos['T'].mean())

        w.append(df_pre['W'] - df_pre['W'].mean())
        w.append(df_ffp['W'] - df_pre['W'].mean())
        w.append(df_pos['W'] - df_pos['W'].mean())

        wt.append((df_pre['T'] - df_pre['T'].mean()) * (df_pre['W'] - df_pre['W'].mean()))
        wt.append((df_ffp['T'] - df_pre['T'].mean()) * (df_ffp['W'] - df_pre['W'].mean()))
        wt.append((df_pos['T'] - df_pos['T'].mean()) * (df_pos['W'] - df_pos['W'].mean()))

        tmean.append([df_pre['T'].mean(), df_ffp['T'].mean(), df_pos['T'].mean()])
        wmean.append([df_pre['W'].mean(), df_ffp['W'].mean(), df_pos['W'].mean()])
        wtmean.append([((df_pre['T'] - df_pre['T'].mean()) * (df_pre['W'] - df_pre['W'].mean())).mean(),
                       ((df_ffp['T'] - df_ffp['T'].mean()) * (df_ffp['W'] - df_ffp['W'].mean())).mean(),
                       ((df_pos['T'] - df_pos['T'].mean()) * (df_pos['W'] - df_pos['W'].mean())).mean()])
    data = [t, w, wt]
    return tmean, wmean, wtmean, data


def get_plot(tmean, wmean, wtmean, data):
    global tower, ti, meantext
    datamean = [tmean, wmean, wtmean]
    color = ['r', 'b', 'lime']
    start = ['15:08:00', '14:55:00', '15:43:00', '14:25:00', '14:55:00']
    ffps = ['15:38:00', '15:25:00', '16:13:00', '14:55:00', '15:25:00']
    ffpe = ['15:53:00', '15:40:00', '16:33:00', '15:15:00', '15:45:00']
    end = ['16:23:00', '16:10:00', '17:03:00', '15:45:00', '16:15:00']

    ylabel = ["t' (\u2103)", "$w' (ms^{-1})$", "$w't' (\u2103ms^{-1})$"]
    yl = -5, -8, -2
    yr = 5, 8, 4
    step = 2, 4, 2

    legend = ['19m AGL', '10m AGL', 'Fire Starts', 'Fire Ends']
    fig, axs = plt.subplots(3, 1, figsize=(10, 6))
    for vi in range(3):
        ax = plt.subplot(3, 1, vi + 1)
        ax.set_ylim(yl[vi], yr[vi])
        plt.yticks(np.arange(yl[vi], yr[vi] + 1, step[vi]), fontsize='large')
        ax.set_xlim(-10, 485)
        ax.set_ylabel(ylabel[vi], fontsize=15)
        plt.subplots_adjust(hspace=0.1, bottom=0.0)
        if vi == 2:
            plt.text(.5, -.22, 'Time', fontsize=15, fontweight='bold', c='k', ha='center', va='center',
                     transform=ax.transAxes)
        for hi in [2, 1]:
            wwpre = data[vi][3 * hi]
            xpre = len(wwpre)

            wwffp = data[vi][3 * hi + 1]
            xffp = len(wwffp)

            wwpos = data[vi][3 * hi + 2]
            xpos = len(wwpos)

            plt.scatter(np.arange(xpre), wwpre, .5, color[hi], alpha=1.)
            plt.scatter(np.arange(xpre, xpre + xffp), wwffp, .5, color[hi], alpha=1., label='_nolegend_')
            plt.scatter(np.arange(xpre + xffp, xpre + xffp + xpos), wwpos, .5, color[hi], alpha=1., label='_nolegend_')

        plt.text(.05, .1, meantext[vi], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
        # plt.text(.12, .1, f'{datamean[vi][1][0]:.3f}', fontsize=12, c=color[1], ha='center', va='center',
        #          transform=ax.transAxes)
        plt.text(.22, .1, f'{datamean[vi][1][0]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.32, .1, f'{datamean[vi][2][0]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)
        # plt.text(.42, .1, f'{datamean[vi][1][1]:.3f}', fontsize=12, c=color[1], ha='center', va='center',
        #          transform=ax.transAxes)
        plt.text(.51, .1, f'{datamean[vi][1][1]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.58, .1, f'{datamean[vi][2][1]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)
        # plt.text(.74, .1, f'{datamean[vi][1][2]:.3f}', fontsize=12, c=color[1], ha='center', va='center',
        #          transform=ax.transAxes)
        plt.text(.84, .1, f'{datamean[vi][1][2]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.94, .1, f'{datamean[vi][2][2]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)

        plt.vlines(xpre, yl[vi], yr[vi], 'magenta', linestyle='--', lw=1)
        plt.vlines(xpre + xffp, yl[vi], yr[vi], 'cyan', linestyle='--', lw=1)
        plt.hlines(0, 0, xpre + xffp + xpos, 'grey', lw=1)
        plt.xticks(np.arange(0, xpre + xffp + xpos, 1000), [])
        if vi == 2:
            plt.text(0, -.08, start[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(.41, -.08, ffps[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(.62, -.08, ffpe[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(1, -.08, end[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
        if vi == 0:
            plt.legend(legend, fontsize=10, markerscale=3., ncol=2, frameon=False, loc='upper right')

    plt.suptitle(f"Control Tower", fontsize=15, fontweight='bold', x=.5, y=.93)
    plt.savefig(f'../../2020_SERDP_Diane/plot/MixedMean/{tower[ti]}/Control_TimeSeries_mixedmean.png', bbox_inches='tight')
    plt.close()

def read_xlsx():
    global tower, ti
    df3 = pd.ExcelFile(f'../../2020_SERDP_Diane/data/PeriodMean/3m/2019_{tower[ti]}_Tower-fft_3mAGL.xlsx')
    df10 = pd.ExcelFile(f'../../2020_SERDP_Diane/data/PeriodMean/10m/2019_{tower[ti]}_Tower-fft_10mAGL.xlsx')
    df20 = pd.ExcelFile(f'../../2020_SERDP_Diane/data/PeriodMean/20m/2019_{tower[ti]}_Tower-fft_20mAGL.xlsx')
    return df3, df10, df20

def get_fluxdata(df3, df10, df20):
    t = []
    w = []
    wt = []
    tmean = []
    wmean = []
    wtmean = []

    df = [df3, df10, df20]
    for hi in range(3):
        df_pre = pd.read_excel(df[hi], 'Pre-FFP')
        df_ffp = pd.read_excel(df[hi], 'FFP')
        df_pos = pd.read_excel(df[hi], 'Post-FFP')

        t.append(df_pre['T (C)'] - df_pre['T (C)'].mean())
        t.append(df_ffp['T (C)'] - df_pre['T (C)'].mean())
        t.append(df_pos['T (C)'] - df_pos['T (C)'].mean())

        w.append(df_pre['W (m/s)'] - df_pre['W (m/s)'].mean())
        w.append(df_ffp['W (m/s)'] - df_pre['W (m/s)'].mean())
        w.append(df_pos['W (m/s)'] - df_pos['W (m/s)'].mean())

        wt.append((df_pre['T (C)'] - df_pre['T (C)'].mean()) * (df_pre['W (m/s)'] - df_pre['W (m/s)'].mean()))
        wt.append((df_ffp['T (C)'] - df_pre['T (C)'].mean()) * (df_ffp['W (m/s)'] - df_pre['W (m/s)'].mean()))
        wt.append((df_pos['T (C)'] - df_pos['T (C)'].mean()) * (df_pos['W (m/s)'] - df_pos['W (m/s)'].mean()))

        tmean.append([df_pre['T (C)'].mean(), df_ffp['T (C)'].mean(), df_pos['T (C)'].mean()])
        wmean.append([df_pre['W (m/s)'].mean(), df_ffp['W (m/s)'].mean(), df_pos['W (m/s)'].mean()])
        wtmean.append(
            [((df_pre['T (C)'] - df_pre['T (C)'].mean()) * (df_pre['W (m/s)'] - df_pre['W (m/s)'].mean())).mean(),
             ((df_ffp['T (C)'] - df_ffp['T (C)'].mean()) * (df_ffp['W (m/s)'] - df_ffp['W (m/s)'].mean())).mean(),
             ((df_pos['T (C)'] - df_pos['T (C)'].mean()) * (df_pos['W (m/s)'] - df_pos['W (m/s)'].mean())).mean()])
    data = [t, w, wt]
    return tmean, wmean, wtmean, data


def get_fluxplot(tmean, wmean, wtmean, data):
    global meantext
    datamean = [tmean, wmean, wtmean]
    color = ['r', 'b', 'lime']
    start = ['15:08:00', '14:55:00', '15:43:00', '14:25:00', '14:55:00']
    ffps = ['15:38:00', '15:25:00', '16:13:00', '14:55:00', '15:25:00']
    ffpe = ['15:53:00', '15:40:00', '16:33:00', '15:15:00', '15:45:00']
    end = ['16:23:00', '16:10:00', '17:03:00', '15:45:00', '16:15:00']

    ylabel = ["t' (\u2103)", "w' (m/s)", "w't' (m\u2103/s)"]
    yl = -5, -8, -100
    yr = 100, 8, 300
    step = 20, 4, 80

    legend = ['19m AGL', '10m AGL', '3m AGL', 'Fire Starts', 'Fire Ends']
    fig, axs = plt.subplots(3, 1, figsize=(10, 6))
    for vi in range(3):
        ax = plt.subplot(3, 1, vi + 1)
        ax.set_ylim(yl[vi], yr[vi])
        plt.yticks(np.arange(yl[vi], yr[vi] + 1, step[vi]), fontsize='large')
        ax.set_xlim(-10, 485)
        ax.set_ylabel(ylabel[vi], fontsize=15)
        plt.subplots_adjust(hspace=0.1, bottom=0.0)
        if vi == 2:
            plt.text(.5, -.22, 'Time', fontsize=15, fontweight='bold', c='k', ha='center', va='center',
                     transform=ax.transAxes)
        for hi in [2, 1, 0]:
            wwpre = data[vi][3 * hi]
            xpre = len(wwpre)

            wwffp = data[vi][3 * hi + 1]
            xffp = len(wwffp)

            wwpos = data[vi][3 * hi + 2]
            xpos = len(wwpos)

            plt.scatter(np.arange(xpre), wwpre, .5, color[hi], alpha=1.)
            plt.scatter(np.arange(xpre, xpre + xffp), wwffp, .5, color[hi], alpha=1., label='_nolegend_')
            plt.scatter(np.arange(xpre + xffp, xpre + xffp + xpos), wwpos, .5, color[hi], alpha=1., label='_nolegend_')

        if vi == 0:
            ht = .9
        else:
            ht = .1
        plt.text(.05, ht, meantext[vi], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
        plt.text(.12, ht, f'{datamean[vi][0][0]:.3f}', fontsize=11, c=color[0], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.22, ht, f'{datamean[vi][1][0]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.32, ht, f'{datamean[vi][2][0]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.435, ht, f'{datamean[vi][0][1]:.3f}', fontsize=11, c=color[0], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.51, ht, f'{datamean[vi][1][1]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.58, ht, f'{datamean[vi][2][1]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.74, ht, f'{datamean[vi][0][2]:.3f}', fontsize=11, c=color[0], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.84, ht, f'{datamean[vi][1][2]:.3f}', fontsize=11, c=color[1], ha='center', va='center',
                 transform=ax.transAxes)
        plt.text(.94, ht, f'{datamean[vi][2][2]:.3f}', fontsize=11, c=color[2], ha='center', va='center',
                 transform=ax.transAxes)

        plt.vlines(xpre, yl[vi], yr[vi], 'magenta', linestyle='--', lw=1)
        plt.vlines(xpre + xffp, yl[vi], yr[vi], 'cyan', linestyle='--', lw=1)
        plt.hlines(0, 0, xpre + xffp + xpos, 'grey', lw=1)
        plt.xticks(np.arange(0, xpre + xffp + xpos, 1000), [])
        if vi == 2:
            plt.text(0, -.08, start[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            #         plt.text(.21, -.08, mid1, fontsize=15, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(.41, -.08, ffps[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(.62, -.08, ffpe[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
            #         plt.text(.82, -.08, mid2, fontsize=15, c='k', ha='center', va='center', transform=ax.transAxes)
            plt.text(1, -.08, end[ti], fontsize=13, c='k', ha='center', va='center', transform=ax.transAxes)
        if vi == 2:
            plt.legend(legend, fontsize=10, markerscale=3., ncol=2, frameon=False)

    plt.suptitle(f"{tower[ti]} Tower", fontsize=15, fontweight='bold', x=.5, y=.93)
    plt.savefig(f'../../2020_SERDP_Diane/plot/MixedMean/{tower[ti]}/{tower[ti]}_TimeSeries_mixedmean.png', bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    tower = ['East', 'Flux', 'North', 'South_Mobile', 'West']
    meantext = ["$\overline{t}$ = ", "$\overline{w}$ = ", "$\overline{w^{'}t^{'}}$= "]

    cstart = [543231, 535431, 564231, 517431, 535431]
    cffps = [561231, 553431, 582231, 535431, 553431]
    cffpe = [570231, 562431, 594231, 547431, 565431]
    cend = [588231, 580431, 612231, 565431, 583431]
    ti = 1

    # df3, df10, df20 = read_table()
    # tmean, wmean, wtmean, data = get_data(df3, df10, df20)
    # get_plot(tmean, wmean, wtmean, data)

    df3, df10, df20 = read_xlsx()
    tmean, wmean, wtmean, data = get_fluxdata(df3, df10, df20)
    get_fluxplot(tmean, wmean, wtmean, data)