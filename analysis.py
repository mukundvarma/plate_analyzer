import numpy as np
import pandas as pd
from itertools import product

VALID_WELLS = [f"{x[0]}{x[1]}" for x in product(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], 1+np.arange(12))]

def collate_results(slopes, intercepts, df):

    results_table = pd.concat([slopes, intercepts], axis=1)
    results_table.columns = ['Rate (RFU/Cycle)', 'Intercept (RFU at t_0)']
    rfu_last = df.iloc[-1, ]
    rfu_last.name = 'End point (RFU at t_last)'
    results_table = pd.merge(results_table, rfu_last, left_index=True, right_index=True, how='right')
    return results_table



def calculate_standard_curve(merged, min_nz, fold_dilutions):
    well_ranks = merged[["Well", "well rank in group"]].value_counts(ascending=True).reset_index(name='count').set_index("Well")["well rank in group"].to_dict()
    means = merged.query('`Data Type` == "Standard Curve"').sort_values("Cycle", ascending=False).head(20).groupby("Well")[
        "RFU (smooth)"
    ].mean().to_frame()
    means['Conc'] = min_nz
    print(means)
    for well in means.index:
        means.loc[well, 'Conc'] = min_nz * (fold_dilutions ** (6 - well_ranks[well]))
    return means

def infer_standard_curve_group(slopes, wells_in_group):
    is_standard_curve = False
    n_low_slopes = 0
    for well in wells_in_group:
        if slopes[well] < 1:
            n_low_slopes += 1

    if n_low_slopes > 2:
        is_standard_curve = True

    return is_standard_curve

def straight_line(x, m, c): # this is your 'straight line' y=f(x)
    return m*x + c

def smooth_concentrations(data, well, delta=5, frac=0.3):
    from statsmodels.nonparametric.smoothers_lowess import lowess
    smoothed = lowess(data[well], data.index, delta=5, frac=0.3)[:, 1]
    return smoothed

def calculate_linear_region(x, y, sensitivity=0.01):
    from kneed import KneeLocator
    
    kl = KneeLocator(x, y, S=sensitivity)
    knee_index_x = np.where(x == kl.knee)[0][0]
    
    if knee_index_x < len(x) // 20:
        knee_index_x = len(x) // 2

    return knee_index_x

def fit_straight_line(x, y):
    from scipy.stats import linregress
    m, c, r_value, p_value, std_err = linregress(x, y)
    y_pred = m * x + c
    return y_pred, x, m, c, r_value

def subtract_background(raw, background_wells=None):
    
    bkg = raw[background_wells].copy()
    subtracted = raw.copy()
    i = 0
    while i < raw.shape[1]:

        subtracted.iloc[:, i] = subtracted.iloc[:, i] - bkg.iloc[:, i // 6] 
        i += 1

    return subtracted
            

def create_merged_df(data, smoothed, fits, slopes):
    raw = pd.melt(data, ignore_index=False)
    smoothed = pd.melt(smoothed, ignore_index=False)
    raw.columns = ['Well', 'raw value']
    smoothed.columns = ['Well', 'smoothed value']
    raw.reset_index(inplace=True)
    smoothed.reset_index(inplace=True)
    raw['id'] = raw['Well'].astype(str) + "_" + raw['Cycle'].astype(str)
    smoothed['id'] = smoothed['Well'].astype(str) + "_" + smoothed['Cycle'].astype(str)
    
    merged = pd.merge(raw, smoothed, left_on='id', right_on='id', how='left')
    merged.index = merged['id']
    merged = merged[['Cycle_x', 'Well_x', 'raw value', 'smoothed value']]
    merged.columns = ['Cycle', 'Well', 'RFU', 'RFU (smooth)']
    merged['Linear Fit'] = np.nan
    merged.reset_index(inplace=True)
    all_fits = {}
    for well, fit in fits.items():
        fit['id'] = [well + "_" +str(s) for s in fit.index]
        fm = fit.set_index('id')['y'].to_dict()
        all_fits = {**all_fits, **fm}
    merged['Linear Fit'] = merged['id'].map(all_fits)
    well_groups = [data.columns[6*i:6*i+6] for i in range(16)]
    merged["Data Type"] = "Reaction"
    merged["Index in group"] = 0
    for group in well_groups:
        if infer_standard_curve_group(slopes, group[:-1]):
            group_idx = merged.query('Well in @group').index
            merged.loc[group_idx, "Data Type"] = "Standard Curve"
        cur_idx = 0
        for well in group:
            cur_idx += 1
            well_idx = merged.query('Well == @well').index
            merged.loc[well_idx, "well rank in group"] = cur_idx

    merged.set_index('id', inplace=True)
    merged['group'] = np.nan
    for i, group in enumerate(well_groups):
        idxs = merged.query('Well in @group').index
        merged.loc[idxs, 'group'] = f"{group[0]} to {group[-1]}"
    return merged



def analyze(raw, background_wells=None):
    assert set(raw.columns) == set(VALID_WELLS)
    if background_wells is None:
        background_wells = raw.columns[5::6]
    
    data = subtract_background(raw, background_wells)
    foreground_wells = data.columns.difference(background_wells)
    foreground_wells = [x for x in data.columns if x in foreground_wells]
    smooth_each = []
    for well in data.columns:
        smooth_each.append(pd.Series(smooth_concentrations(data, well), index=data.index))
    smoothed = pd.concat(smooth_each, axis=1)
    smoothed.index = data.index
    smoothed.columns = data.columns
    
    fits = {}
    slopes = {}
    intercepts = {}
    
    for well in foreground_wells:
        i_linmax = calculate_linear_region(smoothed.index, smoothed[well])
        y_pred, x, m, c, r = fit_straight_line(smoothed.index[:i_linmax], smoothed[well].values[:i_linmax])
        slopes[well] = m
        intercepts[well] = c
        fits[well] = pd.DataFrame(dict(y=y_pred, x=x))
    
    slopes = pd.Series(slopes).loc[foreground_wells]
    intercepts = pd.Series(intercepts).loc[foreground_wells]
    
    merged = create_merged_df(data, smoothed, fits, slopes)

    return merged, slopes, intercepts
