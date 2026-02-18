import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
import sympy as sp
from scipy.optimize import curve_fit
from scipy.optimize import bisect
from collections import defaultdict
from scipy.integrate import quad
from scipy.special import factorial

Dec_G = np.deg2rad(27.12825)
l_NCP = np. deg2rad(122.93192)


D_E = [5800,5550,4500,4250,4150,4100,4100,4100,4100,4100,4100,4100,4000,3900,3800,3750,3650,3600,3550,3500,3450,3350,3250,3150,3050,2950,2800,2650,2500,2300,2100,1900,1750,1500,1250,1050,800,500,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50] # Directional Exposure


def Poission(n,b):
    P=np.exp(-1*b)*np.power(b,n)/np.math.factorial(n)
    return P

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
"""
def zech_integrand(s, n0, b0):
    lam = b0 + s
    return (lam ** n0) * np.exp(-lam) / factorial(n0)
"""
def zech_integrand(s, n0, b0):
    n0 = n0 *43464*81/(118013*225)
    b0 = b0 *43464*81/(118013*225)
    lam = b0 + s
    if lam <= 0:
        return 0.0
    # log-space: lam**n * exp(-lam) / n!  =  exp(n*log(lam) - lam - log(n!))
    log_val = n0 * np.log(lam) - (lam+n0)
    # For extreamly small value put it as zero.
    if log_val < -745:   # np.exp(-745) ~ 5e-324 (lower limit of double)
        return 0.0
    return float(np.exp(log_val))

def zech_upper_limit(n0, b0):
    CL = 0.95
    # Normalize: total integral over [0, ∞)
    Z, _ = quad(zech_integrand, 0, np.inf, args=(n0, b0))

    def integral_diff(sul):
        I, _ = quad(zech_integrand, 0, sul, args=(n0, b0))

        return I / Z - CL

    S_UL = bisect(integral_diff, 0.0, 100.0)
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


df3E1=pd.read_csv('Cut3EeV.csv')
E3L=df3E1['L'].values
E3B=df3E1['B'].values

df3E2=pd.read_csv('test3EeV.csv')
E3bL=df3E2['L'].values
E3bB=df3E2['B'].values
E3bdL=df3E2['dL'].values
E3bdB=df3E2['dB'].values

df3E3=pd.read_csv('JF12_3EeVCut.csv')
E3dL=df3E3['dL'].values
E3dB=df3E3['dB'].values


Dec3E = Declination(E3L,E3B)
bDec3E = Declination(E3bL,E3bB) 
dDec3E = Declination(E3dL,E3dB)
bdDec3E = Declination(E3bdL,E3bdB)



nside=9
npix=hp.nside2npix(nside) #12*nside**2
l_arr,b_arr = get_lb_arrays(nside)
Dec_arr = Declination(l_arr,b_arr)
print("Min Max")
print(np.min(Dec_arr),np.max(Dec_arr))
Dir_Ex = get_DE(Dec_arr)
E3ghpmap=np.zeros(npix)
E3dghpmap=np.zeros(npix)
E3ghpmapb=np.zeros(npix)
E3dghpmapb=np.zeros(npix)

E3gindices=hp.ang2pix(nside,np.pi/2-E3B,E3L)
E3dgindices=hp.ang2pix(nside,np.pi/2-E3dB,E3dL)

E3gindicesb=hp.ang2pix(nside,np.pi/2-E3bB,E3bL)
E3dgindicesb=hp.ang2pix(nside,np.pi/2-E3bdB,E3bdL)

for sig in range(len(E3L)):
    E3ghpmap[E3gindices[sig]]=E3ghpmap[E3gindices[sig]]+1

for back in range(len(E3bL)):
    E3ghpmapb[E3gindicesb[back]]=E3ghpmapb[E3gindicesb[back]]+1

for sig in range(len(E3dL)):
    E3dghpmap[E3dgindices[sig]]=E3dghpmap[E3dgindices[sig]]+1

for back in range(len(E3bdL)):
    E3dghpmapb[E3dgindicesb[back]]=E3dghpmapb[E3dgindicesb[back]]+1




E2gULmap_l=[]
E2dgULmap_l=[]

E3gULmap_l=[]
E3dgULmap_l=[]

"""
for k in range(len(E2dghpmap)):
    if(E2dghpmapb[k]>0):
        b=E2dghpmapb[k]/20
        n=E2dghpmap[k]
        UL=bisect(lambda s: ((s+b)/b)**n - 0.05 * sp.exp(s),0,1000)
        #print(n,b,UL)
        if(b != 0):
            FUL = UL / Dir_Ex[k]
            Dec = np.rad2deg(Dec_arr[k])
            tu = (Dec,FUL)
            E2dgULmap_l.append(tu)
"""
print("NEXT")
            
for k in range(len(E3dghpmap)):
    if(E3dghpmapb[k]>0):
        b=E3dghpmapb[k]/20
        n=E3dghpmap[k]
        UL = zech_upper_limit(n, b)
        #print(n,b,UL)
        if(b != 0 and Dir_Ex[k] > 1200 and UL > 0.000000001):
            FUL = UL / Dir_Ex[k]
            Dec = np.rad2deg(Dec_arr[k])
            tu = (Dec,FUL)
            E3dgULmap_l.append(tu)
print("This is max")
print(np.max(E3dgULmap_l))
Dec3_l, FUL3_l = compute_bin_averages(E3dgULmap_l)

print(Dec3_l[23], FUL3_l)

Um21 = np.mean(FUL3_l[0:24]) * 1000
U15 = np.mean(FUL3_l) * 1000

label = (
    fr"Dec < -21$^\circ$ : {Um21:.2f} $\times 10^{{-3}}$"
    "\n"
    fr"Dec < 15$^\circ$ : {U15:.2f} $\times 10^{{-3}}$"
)

f=open("EvsR_3EeV.csv",'w')

head = "Dec, FUL \n"
f.write(head)
for d in  range(len(Dec3_l)):
    data = "%lf, %lf \n" %(Dec3_l[d], FUL3_l[d])
    f.write(data)
f.close()
plt.figure(figsize=(12, 8))
plt.title('Mean Flux Upper Limit, E > 3 EeV, nside = 9', fontsize=30)
plt.step(Dec3_l, FUL3_l, where='mid', color='red', label=label, linewidth=1.5)

print(np.mean(FUL3_l[0:24]), np.mean(FUL3_l))

# Set axis labels
plt.xlabel(r'Declination [°]', fontsize=25)
plt.ylabel(r'Mean Flux Upper Limit [$10^{-3} \, \mathrm{km}^{-2} \, \mathrm{yr}^{-1}$]', fontsize=25)

# Format y-axis in units of 10^-3
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x * 1e3:.0f}'))

# Add legend
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)

# Show the plot
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()
