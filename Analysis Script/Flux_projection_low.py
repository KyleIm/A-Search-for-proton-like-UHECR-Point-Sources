import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from scipy.integrate import quad


D_E = [5800,5750,4600,4250,4150,
         4100,4100,4100,4100,4000,
         4000,3950,3900,3800,3700,
         3600,3500,3450,3400,3300,
         3250,3150,3100,3050,2950,
         2850,2700,2600,2450,2300,
         2100,1900,1750,1500,1250,
         1050,800,500,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50] 

Dec_G = np.deg2rad(27.12825)
l_NCP = np. deg2rad(122.93192)

def Declination(l,b):
    l=np.deg2rad(l)
    b=np.deg2rad(b)
    Dec = np.arcsin(np.sin(b)*np.sin(Dec_G) + np.cos(b)*np.cos(Dec_G)*np.cos(l_NCP- l))
    return Dec

def get_DE(Dec):
    declination_rad = np.asarray(Dec)
    declination_deg = np.degrees(declination_rad)
    indices = np.clip(((declination_deg + 90) // 3).astype(int), 0, len(D_E) - 1)
    return np.array(D_E)[indices]


def get_lb_arrays(nside):
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    b_array = np.degrees(0.5 * np.pi - theta)
    l_array = np.degrees(phi)    
    return l_array, b_array

def zech_integrand(s, n0, b0):
    n0 = n0 *43464*36/(175952*225)
    b0 = b0 *43464*36/(175952*225)
    lam = b0 + s
    if lam <= 0:
        return 0.0
    log_val = n0 * np.log(lam) - (lam+n0)
    if log_val < -745:
        return 0.0
    return float(np.exp(log_val))


def zech_upper_limit(n0, b0):
    CL = 0.95
    Z, _ = quad(zech_integrand, 0, np.inf, args=(n0, b0))

    def integral_diff(sul):
        I, _ = quad(zech_integrand, 0, sul, args=(n0, b0))
        return I / Z - CL

    from scipy.optimize import bisect
    S_UL = bisect(integral_diff, 0.0, 100.0)
    return S_UL


def format_power_text(value):
    if value <= 0 or not np.isfinite(value):
        return ""
    exponent = int(np.floor(np.log10(value)))
    coeff = value / (10 ** exponent)
    return rf"${coeff:.2f}\times10^{{{exponent}}}$"


def format_major_tick(value, _pos):
    if value <= 0 or not np.isfinite(value):
        return ""
    exponent = int(np.floor(np.log10(value)))
    coeff = value / (10 ** exponent)
    rounded_coeff = int(np.round(coeff))
    if np.isclose(coeff, rounded_coeff, atol=1e-10, rtol=0.0):
        coeff_str = f"{rounded_coeff}"
    else:
        coeff_str = f"{coeff:.1f}".rstrip("0").rstrip(".")
    return rf"${coeff_str}\times10^{{{exponent}}}$"

df1=pd.read_csv('JF12_2to3EeVCut.csv')
dL=df1['dL'].values
dB=df1['dB'].values

df2=pd.read_csv('test2to3EeV.csv')
bdL=df2['dL'].values
bdB=df2['dB'].values


nside=6
npix=hp.nside2npix(nside)
ghpmap=np.zeros(npix)
dghpmap=np.zeros(npix)
dgindices=hp.ang2pix(nside,np.pi/2-dB,dL)

l_arr,b_arr = get_lb_arrays(nside)
Dec_arr = Declination(l_arr,b_arr)
Dir_Ex = get_DE(Dec_arr)
dghpmapb=np.zeros(npix)
dgindicesb=hp.ang2pix(nside,np.pi/2-bdB,bdL)

for sig in range(len(dL)):
    dghpmap[dgindices[sig]]=dghpmap[dgindices[sig]]+1

for back in range(len(bdL)):
    dghpmapb[dgindicesb[back]]=dghpmapb[dgindicesb[back]]+1

dgULmap_l=[]

for k in range(len(dghpmap)):
    b=dghpmapb[k]/20
    n=dghpmap[k]
    if(b==0):
        FUL=np.nan
    else:
        if(Dir_Ex[k] > 1200):
            UL = zech_upper_limit(n, b)
            FUL =UL/Dir_Ex[k]
            if(UL < 0.000000001):
                print(UL,FUL,Dir_Ex[k])
                FUL = np.nan
        else:
            FUL = 100
    dgULmap_l.append(FUL)

hpmap = np.array(dgULmap_l)
masked_hpmap = np.where(hpmap < 99, hpmap, np.nan)
positive_vals = masked_hpmap[np.isfinite(masked_hpmap) & (masked_hpmap > 0)]
if positive_vals.size == 0:
    raise ValueError("No positive finite values available for log-scale projection.")
vmin = float(np.nanmin(positive_vals))
vmax = float(np.nanmax(positive_vals))
if np.isclose(vmin, vmax, rtol=0.0, atol=1e-15):
    vmax = vmin * 10
        
projview(
    masked_hpmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,    
    xlabel="L",
    ylabel="B",
    flip="astro",
    projection_type="mollweide",
    norm=LogNorm(vmin=vmin, vmax=vmax)
)

ax = plt.gca()
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(20)

cbar = plt.gcf().axes[-1]
base = np.round(vmin / (10 ** np.floor(np.log10(vmin)))) * (10 ** np.floor(np.log10(vmin)))
if not np.isfinite(base) or base <= 0:
    base = vmin
major_candidates = base * np.arange(1, int(np.floor(vmax / base)) + 1)
major_candidates = major_candidates[(major_candidates > vmin) & (major_candidates < vmax)]
major_candidates = major_candidates[~np.isclose(major_candidates, 1e-3, atol=1e-12, rtol=0.0)]

if major_candidates.size:
    log_gap_min = np.abs(np.log10(major_candidates) - np.log10(vmin))
    log_gap_max = np.abs(np.log10(vmax) - np.log10(major_candidates))
    major_ticks = major_candidates[(log_gap_min > 0.08) & (log_gap_max > 0.08)]
else:
    major_ticks = major_candidates

is_vertical = cbar.get_position().height >= cbar.get_position().width
if is_vertical:
    cbar.set_yticks(major_ticks)
    cbar.yaxis.set_major_formatter(FuncFormatter(format_major_tick))
    cbar.tick_params(axis="y", labelsize=14)
    cbar.yaxis.offsetText.set_visible(False)
    cbar.text(0.5, -0.08, format_power_text(vmin), transform=cbar.transAxes, ha="center", va="top", fontsize=20)
    cbar.text(0.5, 1.08, format_power_text(vmax), transform=cbar.transAxes, ha="center", va="bottom", fontsize=20)
else:
    cbar.set_xticks(major_ticks)
    cbar.xaxis.set_major_formatter(FuncFormatter(format_major_tick))
    cbar.tick_params(axis="x", labelsize=14)
    cbar.xaxis.offsetText.set_visible(False)
    cbar.text(-0.02, -0.55, format_power_text(vmin), transform=cbar.transAxes, ha="right", va="center", fontsize=20)
    cbar.text(1.02, -0.55, format_power_text(vmax), transform=cbar.transAxes, ha="left", va="center", fontsize=20)
cbar.minorticks_off()

plt.title(r"Flux Upper Limit projection, 2EeV  < E $\leq$ 3EeV, $N_{side}$ = 6",fontsize=30)
plt.xlabel("L",fontsize=25)
plt.ylabel("B",fontsize=25)
plt.show()
