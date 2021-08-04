import pandas as pd


def get_data(ht, prd):
    xls = pd.ExcelFile(f'./data/PeriodMean/{ht}/2019_East_Tower-fft_{ht}AGL.xlsx')
    df_in = pd.read_excel(xls, prd)
    df_out = pd.DataFrame({})
    df_out['Time'] = df_in.iloc[:, 0].dropna()
    df_out[f'W_{ht}'] = df_in['W (m/s)'].dropna()
    df_out[f'T_{ht}'] = df_in['T (C)'].dropna()
    return df_out


def merge_btw_ht(prd):
    df3 = get_data('3m', prd)
    df10 = get_data('10m', prd)
    df20 = get_data('20m', prd)
    df_temp = df20.merge(df10, on=['Time'], how='outer')[['Time', 'W_20m', 'T_20m', 'W_10m', 'T_10m']]
    df = df_temp.merge(df3, on=['Time'], how='outer')[['Time', 'W_20m', 'T_20m', 'W_10m', 'T_10m', 'W_3m', 'T_3m']]
    return df


if __name__ == "__main__":
    prd = ['Pre-FFP', 'FFP', 'Post-FFP']
    for p in prd:
        df = merge_btw_ht(p)
        df.to_csv(f'./data/East_{p}.csv', index=False)

    # df[df.isnull().any(axis=1)] ## locate nan values
