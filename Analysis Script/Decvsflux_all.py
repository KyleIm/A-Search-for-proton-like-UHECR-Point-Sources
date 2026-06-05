import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from scipy.optimize import bisect
from scipy.integrate import quad
from matplotlib.ticker import FuncFormatter
from scipy.special import gammaincc

Dec_G = np.deg2rad(27.12825)
l_NCP = np. deg2rad(122.93192)


D_E = [5800,5550,4500,4250,4150,4100,4100,4100,4100,4100,4100,4100,4000,3900,3800,3750,3650,3600,3550,3500,3450,3350,3250,3150,3050,2950,2800,2650,2500,2300,2100,1900,1750,1500,1250,1050,800,500,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50] # Directional Exposure

D_E_low = [5800,5750,4600,4250,4150,
           4100,4100,4100,4100,4000,
           4000,3950,3900,3800,3700,
           3600,3500,3450,3400,3300,
           3250,3150,3100,3050,2950,
           2850,2700,2600,2450,2300,
           2100,1900,1750,1500,1250,
           1050,800,500,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50] # Directional Exposure for 2to3EeV

def Declination(l,b):
    l=np.deg2rad(l)
    b=np.deg2rad(b)
    Dec = np.arcsin(np.sin(b)*np.sin(Dec_G) + np.cos(b)*np.cos(Dec_G)*np.cos(l_NCP- l))
    return Dec

def get_DE(Dec, exposure_table):
    declination_rad = np.asarray(Dec)
    declination_deg = np.degrees(declination_rad)
    indices = np.clip(((declination_deg + 90) // 3).astype(int), 0, len(exposure_table) - 1)
    return np.array(exposure_table)[indices]


def get_lb_arrays(nside):
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    b_array = np.degrees(0.5 * np.pi - theta)
    l_array = np.degrees(phi)    
    return l_array, b_array
"""
def zech_integrand(s, n0, b0):
    lam = b0 + s
    return (lam ** n0) * np.exp(-lam) / factorial(n0)
"""
def zech_integrand(s, n0, b0, norm_factor=1.0):
    n0 = n0 * norm_factor
    b0 = b0 * norm_factor
    lam = b0 + s
    if lam <= 0:
        return 0.0
    # log-space: lam**n * exp(-lam) / n!  =  exp(n*log(lam) - lam - log(n!))
    log_val = n0 * np.log(lam) - (lam+n0)
    # For extreamly small value put it as zero.
    if log_val < -745:   # np.exp(-745) ~ 5e-324 (lower limit of double)
        return 0.0
    return float(np.exp(log_val))

def zech_upper_limit(n0, b0, norm_factor=1.0):
    CL = 0.95
    n_eff = n0 * norm_factor
    b_eff = b0 * norm_factor
    alpha = n_eff + 1.0

    def integral_diff(sul):
        # I/Z = [∫_b^(b+s) t^n e^-t dt] / [∫_b^∞ t^n e^-t dt]
        #     = [Q(alpha,b) - Q(alpha,b+s)] / Q(alpha,b), Q=gammaincc
        q_b = gammaincc(alpha, b_eff)
        q_bs = gammaincc(alpha, b_eff + sul)
        if q_b <= 0 or not np.isfinite(q_b) or not np.isfinite(q_bs):
            return np.nan
        return (q_b - q_bs) / q_b - CL

    low = 0.0
    high = 100.0
    f_low = integral_diff(low)
    f_high = integral_diff(high)

    if not np.isfinite(f_low):
        return np.nan

    # Expand upper bracket until the root is enclosed.
    while (not np.isfinite(f_high) or f_high <= 0.0) and high < 1e7:
        high *= 2.0
        f_high = integral_diff(high)

    if not np.isfinite(f_high) or f_high <= 0.0:
        return np.nan

    S_UL = bisect(integral_diff, low, high)
    return S_UL



def compute_bin_averages(data):
    bin_size=3
    dec_min=-90
    dec_max=15
    """
    Compute the average ful values for specified declination bins.
    
    Parameters:
        data (list of tuples): Input list in the form [(deg1, ful1), (deg2, ful2), ..., (degN, fulN)].
        bin_size (int): Size of each declination bin in degrees.
        dec_min (int): Minimum declination value.
        dec_max (int): Maximum declination value.
    
    Returns:
        list: Average ful values for each bin.
    """
    # Define bins and bin centers
    bins = np.arange(dec_min, dec_max + bin_size, bin_size)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    # Initialize list for average ful values
    ful_averages = []
    
    # Loop through bins and calculate the average ful for each bin
    for bin_left, bin_right in zip(bins[:-1], bins[1:]):
        # Extract ful values for the current bin
        ful_values = [ful for dec, ful in data if bin_left < dec <= bin_right]
        # Compute the average ful or set to np.nan if no values exist
        ful_averages.append(np.mean(ful_values) if ful_values else np.nan)
    
    return bin_centers, ful_averages


def build_dec_flux_curve(signal_csv, background_csv, nside, exposure_table, norm_factor=1.0):
    df1 = pd.read_csv(signal_csv)
    dL = df1['dL'].values
    dB = df1['dB'].values

    df2 = pd.read_csv(background_csv)
    bdL = df2['dL'].values
    bdB = df2['dB'].values

    npix = hp.nside2npix(nside)
    l_arr, b_arr = get_lb_arrays(nside)
    dec_arr = Declination(l_arr, b_arr)
    dir_ex = get_DE(dec_arr, exposure_table)

    dghpmap = np.zeros(npix)
    dghpmapb = np.zeros(npix)
    dgindices = hp.ang2pix(nside, np.pi / 2 - dB, dL)
    dgindicesb = hp.ang2pix(nside, np.pi / 2 - bdB, bdL)

    for idx in dgindices:
        dghpmap[idx] += 1

    for idx in dgindicesb:
        dghpmapb[idx] += 1

    dec_flux_pairs = []
    um21_vals = []
    u15_vals = []

    for k in range(len(dghpmap)):
        if dghpmapb[k] <= 0:
            continue
        b = dghpmapb[k] / 20
        n = dghpmap[k]
        ul = zech_upper_limit(n, b, norm_factor)
        if b != 0 and dir_ex[k] >= 1250 and ul > 1e-9:
            ful = ul / dir_ex[k]
            dec = np.rad2deg(dec_arr[k])
            dec_flux_pairs.append((dec, ful))
            u15_vals.append(ful)
            if dir_ex[k] >= 3150:
                um21_vals.append(ful)

    dec_l, ful_l = compute_bin_averages(dec_flux_pairs)
    um21 = np.mean(um21_vals) * 1000 if um21_vals else np.nan
    u15 = np.mean(u15_vals) * 1000 if u15_vals else np.nan
    return dec_l, ful_l, um21, u15


datasets = [
    {
        "signal_csv": "JF12_5EeVCut.csv",
        "background_csv": "test5EeV.csv",
        "nside": 15,
        "color": "red",
        "title": r"E > 5 EeV",
        "exposure_table": D_E,
        "norm_factor": 1.0,
    },
    {
        "signal_csv": "JF12_3EeVCut.csv",
        "background_csv": "test3EeV.csv",
        "nside": 9,
        "color": "blue",
        "title": r"E > 3 EeV",
        "exposure_table": D_E,
        "norm_factor": (43464 * 81) / (118013 * 225),
    },
    {
        "signal_csv": "JF12_2to3EeVCut.csv",
        "background_csv": "test2to3EeV.csv",
        "nside": 6,
        "color": "green",
        "title": r"2 EeV < E $\leq$ 3 EeV",
        "exposure_table": D_E_low,
        "norm_factor": (43464 * 81) / (175952 * 225),
    },
]

plt.figure(figsize=(12, 8))
plt.title(r'Mean Flux Upper Limit vs Declination', fontsize=30)

for cfg in datasets:
    dec_l, ful_l, um21, u15 = build_dec_flux_curve(
        cfg["signal_csv"],
        cfg["background_csv"],
        cfg["nside"],
        cfg["exposure_table"],
        cfg["norm_factor"],
    )
    label = (
        fr"{cfg['title']} ($N_{{side}}$={cfg['nside']}) | "
        fr"Dec < -21$^\circ$: {um21:.2f} $\times 10^{{-3}}$, "
        fr"Dec < 15$^\circ$: {u15:.2f} $\times 10^{{-3}}$"
    )
    plt.step(dec_l, ful_l, where='mid', color=cfg["color"], label=label, linewidth=1.5)

# Set axis labels
plt.xlabel(r'Declination [°]', fontsize=25)
plt.ylabel(r'Mean Flux Upper Limit [$10^{-3} \, \mathrm{km}^{-2} \, \mathrm{yr}^{-1}$]', fontsize=25)

# Format y-axis in units of 10^-3
plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x * 1e3:.0f}'))

# Add legend
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)

# Show the plot
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()
