import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from scipy.optimize import curve_fit

df1=pd.read_csv('Cut5EeV.csv')
ID=df1['AugerID'].values
L=df1['L'].values
B=df1['B'].values

df2=pd.read_csv('test5EeV.csv')
Trial=df2['Trial'].values
bL=df2['L'].values
bB=df2['B'].values
bdL=df2['dL'].values
bdB=df2['dB'].values

df3=pd.read_csv('JF12_5EeVCut.csv')
ID3=df3['AugerID'].values
dL=df3['dL'].values
dB=df3['dB'].values

df4=pd.read_csv('Realpsf.csv')
uc_map = df4['uc_psf_map'].values
jf12_map = df4['psf_map'].values

df5=pd.read_csv('Background_psf.csv')
ucback = df5['ucback_psf_map'].values
back=df5['back_psf_map'].values


def Gaussian(x,A,mean,std):
    fit_line=A*np.exp(-0.5*((x-mean)/(std*np.sqrt(std)))**2)
    return fit_line

def LiMa(n,b,alpha):
    if(b==n):
        lima=np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    elif(n==0):
        lima=1234567
    else:
        lima=(n-b)/abs(n-b)*np.sqrt(2*(n*np.log((n+alpha*n)/(b+alpha*n))+b/alpha*np.log((b+alpha*b)/(b+alpha*n))))
    return lima

def Angdistance(l1,l2,b1,b2): #Analytic formula for distance on sphere
    D=np.arccos(np.cos(b1)*np.cos(b2) + np.sin(b1)*np.sin(b2)*np.cos(l2-l1))
    return D

def PSF(D, sigma):
    psf = np.exp(-1*D**2/(2*sigma**2))/(2*np.pi*sigma**2)
    return psf

ucpsf_map = uc_map/np.sum(uc_map) * (43463) * (125*125/225)
psf_map = jf12_map/np.sum(jf12_map) * (43463)  * (125*125/225)
ucbackpsf_map = ucback/np.sum(ucback) * (43463*20)  * (125*125/225)
backpsf_map = back/np.sum(back) * (43463*20)  * (125*125/225)

print("START")

    
nside=128
sigma = np.deg2rad(2.16)
npix=hp.nside2npix(nside)

ghpmap=ucpsf_map
dghpmap=psf_map

ghpmapb=ucbackpsf_map
dghpmapb=backpsf_map

glimamap_l=[]
dglimamap_l=[]

print("Pass")

for k in range(len(ghpmap)):
    if(ghpmapb[k]<0.01):
        pass
    else:
        lima=LiMa(ghpmap[k],ghpmapb[k]/20,1/20)
        if(lima != 1234567):
            glimamap_l.append(lima)

for k in range(len(dghpmap)):
    if(dghpmapb[k]<0.01):
        pass
    else:
        lima=LiMa(dghpmap[k],dghpmapb[k]/20,1/20)
        if(lima != 1234567):
            dglimamap_l.append(lima)

print(np.min(dglimamap_l),np.min(glimamap_l),np.max(dglimamap_l),np.max(glimamap_l))
print(len(glimamap_l))

H_l,bin_edges=np.histogram(glimamap_l,bins=30,range=(-2.5,2.5))


sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1


x_l=np.linspace(-2.5,2.5,30)
ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)


gmean=np.mean(glimamap_l)
gstd=np.std(glimamap_l)
gxval=np.linspace(-3.75,3.75,150)
gyval=[]
for i in gxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    #y=Gaussian(i,gmean,gstd)
    gyval.append(y)

print(ppot[0],ppot[1],ppot[2])

binsize=int(np.max(dglimamap_l))
H_l,bin_edges=np.histogram(dglimamap_l,bins=30,range=(-2,2))
sigma_l=np.sqrt(H_l)

for j in range(len(sigma_l)):
    if(sigma_l[j]==0):
        sigma_l[j]=1

dgmean=np.mean(dglimamap_l)
dgstd=np.std(dglimamap_l)
x_l=np.linspace(-2,2,30)
ppot,pcov=curve_fit(Gaussian,x_l,H_l,p0=[150,-0.03,0.844],sigma=sigma_l)
dgxval=np.linspace(-2.5,2.5,150)
dgyval=[]

print(ppot[0],ppot[1],ppot[2])

for i in dgxval:
    y=Gaussian(i,ppot[0],ppot[1],ppot[2])
    dgyval.append(y)

H_raw, edges = np.histogram(glimamap_l, bins=30, range=(-0.25, 0.25))
centers = 0.5 * (edges[:-1] + edges[1:])
mask = H_raw > 0
p0_raw = [H_raw.max(), np.mean(glimamap_l), np.std(glimamap_l)]
pp_raw, pcov_raw = curve_fit(
    Gaussian,
    centers[mask], H_raw[mask],
    p0=p0_raw,
    sigma=np.sqrt(H_raw[mask]),
    absolute_sigma=False,
    maxfev=10000
)


H_jf, edges_jf = np.histogram(dglimamap_l, bins=30, range=(-0.25, 0.25))
centers_jf = 0.5 * (edges_jf[:-1] + edges_jf[1:])
mask_jf = H_jf > 0
p0_jf = [H_jf.max(), np.mean(dglimamap_l), np.std(dglimamap_l)]
pp_jf, pcov_jf = curve_fit(
    Gaussian,
    centers_jf[mask_jf], H_jf[mask_jf],
    p0=p0_jf,
    sigma=np.sqrt(H_jf[mask_jf]),
    absolute_sigma=False,
    maxfev=10000
)

    
plt.figure()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


xfit = np.linspace(-2.5, 2.5, 400)
#plt.plot(xfit, Gaussian(xfit, *pp_raw), lw=3, color='blue')
#plt.plot(xfit, Gaussian(xfit, *pp_jf), lw=2, color='red')


plt.xlabel("Li-Ma Significance",fontsize = 25)
plt.hist(glimamap_l,bins=30,range=(-2.5,2.5),facecolor='None',lw=3,edgecolor='blue',label='RAW : %d'%len(glimamap_l))
plt.hist(dglimamap_l,bins=30,range=(-2.5,2.5),facecolor='None',lw=2,edgecolor='red',label='JF12 : %d'%len(dglimamap_l))
plt.legend(fontsize=20)
#plt.text(-2.7,ppot[0]*1.2,"A : %0.2f\n$\mu$ : %0.4f \n$\sigma$ : %0.3f"%(ppot[0],ppot[1],ppot[2]),fontsize=25,bbox=dict(boxstyle="square",facecolor='None',edgecolor='black'))
plt.title("Li-Ma : RAW vs JF12, bin=2700, alpha=1/20, E > 5EeV - log (with PSF)", fontsize=30)
plt.yscale('log')
plt.show()



#I need to plot Gaussian Fitting plot
