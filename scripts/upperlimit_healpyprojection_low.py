import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview
from scipy.optimize import bisect
import sympy as sp
from matplotlib.colors import LogNorm
from scipy.integrate import quad
from scipy.special import factorial
from scipy.stats import poisson
"""
D_E = [5800,5550,4500,4250,4150,4100,4100,4100,4100,4100,4100,4100,4000,3900,3800,3750,3650,3600,3550,3500,3450,3350,3250,3150,3050,2950,2800,2650,2500,2300,2100,1900,1750,1500,1250,1050,800,500,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50] # Directional Exposure
"""
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

    from scipy.optimize import bisect
    S_UL = bisect(integral_diff, 0.0, 100.0)
    return S_UL

"""
df1=pd.read_csv('Cut3EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=df1['B'].values

df2=pd.read_csv('test3EeV.csv')
Trial=df2['Trial'].values
bL=df2['L'].values
bB=df2['B'].values
bdL=df2['dL'].values
bdB=df2['dB'].values

df3=pd.read_csv('JF12_3EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=df3['dB'].values
"""
df1=pd.read_csv('Cut2to3EeV.csv')
L=df1['L'].values
B=df1['B'].values

df2=pd.read_csv('test2to3EeV.csv')
bL=df2['L'].values
bB=df2['B'].values
bdL=df2['dL'].values
bdB=df2['dB'].values

df3=pd.read_csv('JF12_2to3EeVCut.csv')
dL=df3['dL'].values
dB=df3['dB'].values

nside=6
npix=hp.nside2npix(nside) #12*nside**2
ghpmap=np.zeros(npix)
dghpmap=np.zeros(npix)
gindices=hp.ang2pix(nside,np.pi/2-B,L)
dgindices=hp.ang2pix(nside,np.pi/2-dB,dL)

l_arr,b_arr = get_lb_arrays(nside)
Dec_arr = Declination(l_arr,b_arr)
Dir_Ex = get_DE(Dec_arr)
ghpmapb=np.zeros(npix)
dghpmapb=np.zeros(npix)
gindicesb=hp.ang2pix(nside,np.pi/2-bB,bL)
dgindicesb=hp.ang2pix(nside,np.pi/2-bdB,bdL)

print("PL")

for sig in range(len(L)):
    ghpmap[gindices[sig]]=ghpmap[gindices[sig]]+1

for back in range(len(bL)):
    ghpmapb[gindicesb[back]]=ghpmapb[gindicesb[back]]+1

for sig in range(len(dL)):
    dghpmap[dgindices[sig]]=dghpmap[dgindices[sig]]+1

for back in range(len(bdL)):
    dghpmapb[dgindicesb[back]]=dghpmapb[dgindicesb[back]]+1

glimamap_l=[]
dgULmap_l=[]

s= sp.symbols('s')

print("StartLOOP")

print(dghpmap)

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
                UL = np.nan
                print("Wrong bin : %d"%k)
        else:
            FUL = 100
    dgULmap_l.append(FUL) #Temperorly change for Cerca

print(len(dgULmap_l))
#18 23 16

hpmap = np.array(dgULmap_l)
#masked_hpmap = np.where(hpmap <0.05991464547104641, hpmap, np.nan)
masked_hpmap = np.where(hpmap < 99, hpmap, np.nan)
log_hpmap = np.log10(hpmap + 1e-10)
        
projview(
    masked_hpmap, #Tempeorly change for Cerca
    coord=["G"],
    graticule=True,
    graticule_labels=True,    
    xlabel="L",
    ylabel="B",
    flip="astro",
    projection_type="mollweide",
    norm=LogNorm(vmin=2, vmax=615)
)

#print(masked_hpmap[1850], masked_hpmap[2150], masked_hpmap[2200])

ax = plt.gca()
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(20)

cbar = plt.gcf().axes[-1]
cbar.tick_params(labelsize=20)



#print(hpmap,np.max(hpmap))
plt.title(r"Log-Scale Flux Upper Limit projection : 432 Bins, 2EeV <  E $\leq$ 3EeV",fontsize=30)
#plt.title("Log-Scale Flux Upper Limit projection (with Extra Events) : 2700 Bins,  E>5EeV",fontsize=30)
plt.xlabel("L",fontsize=25)
plt.ylabel("B",fontsize=25)
#plt.figure()
#plt.subplot(111,projection='aitoff')
#plt.hexbin(dl,db)
#plt.hist(hpmap,bins=14)
#plt.yscale('log')
plt.show()
