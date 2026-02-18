import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview
from scipy.optimize import curve_fit
import time as ti

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=np.pi/2 - df1['B'].values

df2=pd.read_csv('test5EeV.csv')
bL=df2['L'].values
bB=np.pi/2 - df2['B'].values
dbL=df2['dL'].values
dbB=np.pi/2 - df2['dB'].values


df3=pd.read_csv('JF12_5EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=np.pi/2 - df3['dB'].values

def Gaussian(x,A,mean,std):
    fit_line=A*np.exp(-0.5*((x-mean)/(std*np.sqrt(std)))**2)
    return fit_line

def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
        print("Something went wrong")
    elif(n==0):
        lima=(n-b)/abs(n-b)*np.sqrt(2*(b/alpha*np.log((b+alpha*b)/(b+alpha*n)))) #Not making inside ln tobe 0)
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    if(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))<0):
        print("Imaginary")
    return lima


def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)
    return psf


nside=15
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside)


start = ti.time()


theta, phi = hp.pix2ang(nside, np.arange(npix))

org_psf_map = np.zeros(npix, dtype=float)
back_psf_map = np.zeros(npix, dtype=float)

ucorg_psf_map = np.zeros(npix, dtype=float)
ucback_psf_map = np.zeros(npix, dtype=float)


print(ti.time() - start)


print("This will take time")
for i in range(len(dL)):
    D = Angdistance(dL[i], phi, dB[i], theta)
    org_psf_map += PSF(D, sigma)
    UD = Angdistance(L[i], phi, B[i], theta)
    ucorg_psf_map += PSF(UD, sigma)

for i in range(len(dbL)):
    D = Angdistance(dbL[i], phi, dbB[i], theta)
    back_psf_map += PSF(D, sigma)
    UD = Angdistance(bL[i], phi, bB[i], theta)
    ucback_psf_map += PSF(UD, sigma)

org_psf_map = org_psf_map/np.sum(org_psf_map) * 43464
ucorg_psf_map = org_psf_map/np.sum(ucorg_psf_map) * 43464
print(ti.time()-start)

#Here we should calculate Li-Ma

uclima_map = []
lima_map = []

for j in range(len(ucorg_psf_map)):
    if(back_psf_map[j]<1e-2):
        pass
        #uclima_map.append(np.nan)
    else:
        lima=LiMa(ucorg_psf_map[j],ucback_psf_map[j]/20,1/20)
        uclima_map.append(lima)

for k in range(len(org_psf_map)):
    if(back_psf_map[k]<1e-2):
        pass
        #lima_map.append(np.nan)
    else:
        lima=LiMa(org_psf_map[k],back_psf_map[k]/20,1/20)
        lima_map.append(lima)
print(ti.time()-start)
H_l,bin_edges=np.histogram(uclima_map,bins=30,range=(-3.75,3.75))

"""
sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1


x_l=np.linspace(-3.75,3.75,30)
ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)


gmean=np.mean(uclima_map)
gstd=np.std(uclima_map)
gxval=np.linspace(-3.75,3.75,150)
gyval=[]
for i in gxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    #y=Gaussian(i,gmean,gstd)
    gyval.append(y)

print(ppot[0],ppot[1],ppot[2])

binsize=int(np.max(lima_map))
H_l,bin_edges=np.histogram(lima_map,bins=30,range=(-3.525,3.755))
sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1

dgmean=np.mean(lima_map)
dgstd=np.std(lima_map)
x_l=np.linspace(-3.75,3.75,30)
#ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)
dgxval=np.linspace(-3.75,3.75,150)
dgyval=[]

#print(ppot[0],ppot[1],ppot[2])


for i in dgxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    dgyval.append(y)



H_raw, edges = np.histogram(uclima_map, bins=30, range=(-3.75, 3.75))
centers = 0.5 * (edges[:-1] + edges[1:])
mask = H_raw > 0
p0_raw = [H_raw.max(), np.mean(uclima_map), np.std(uclima_map)]
pp_raw, pcov_raw = curve_fit(
    Gaussian,
    centers[mask], H_raw[mask],
    p0=p0_raw,
    sigma=np.sqrt(H_raw[mask]),
    absolute_sigma=False,
    maxfev=10000
)


H_jf, edges_jf = np.histogram(lima_map, bins=30, range=(-3.75, 3.75))
centers_jf = 0.5 * (edges_jf[:-1] + edges_jf[1:])
mask_jf = H_jf > 0
p0_jf = [H_jf.max(), np.mean(lima_map), np.std(lima_map)]
pp_jf, pcov_jf = curve_fit(
    Gaussian,
    centers_jf[mask_jf], H_jf[mask_jf],
    p0=p0_jf,
    sigma=np.sqrt(H_jf[mask_jf]),
    absolute_sigma=False,
    maxfev=10000
)
"""

    
plt.figure()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#xfit = np.linspace(-3.75, 3.75, 400)
#plt.plot(xfit, Gaussian(xfit, *pp_raw), lw=3, color='blue')
#plt.plot(xfit, Gaussian(xfit, *pp_jf), lw=2, color='red')

plt.xlabel("Li-Ma Significance",fontsize = 25)
plt.hist(uclima_map,bins=30,range=(-3.75,3.75),facecolor='None',lw=3,edgecolor='blue',label='RAW : %d'%len(uclima_map))
plt.hist(lima_map,bins=30,range=(-3.75,3.75),facecolor='None',lw=2,edgecolor='red',label='JF12 : %d'%len(lima_map))
plt.legend(fontsize=20)
#plt.text(-3.9,ppot[0]*0.35,"A : %0.2f\n$\mu$ : %0.4f \n$\sigma$ : %0.3f"%(ppot[0],ppot[1],ppot[2]),fontsize=25,bbox=dict(boxstyle="square",facecolor='None',edgecolor='black'))
plt.title("Li-Ma : RAW vs JF12(random), bin=2700, alpha=1/20, E > 5EeV - log", fontsize=30)
plt.yscale('log')
plt.show()
