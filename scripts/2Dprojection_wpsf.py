import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df=pd.read_csv('Realpsf.csv')
psfmap = df['uc_psf_map'].values

nside=128
hpmap=psfmap/np.sum(psfmap) * 43463 * 128 * 128 / 225


projview(
    hpmap,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    #title="Binning with healpix",
    xlabel="l",
    ylabel="b",
    flip="astro",
    projection_type="mollweide",
)

ax = plt.gca()
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(20)  

ax.set_xlabel("L", fontsize=25)
ax.set_ylabel("B", fontsize=25)
    
cbar = plt.gcf().axes[-1]
cbar.tick_params(labelsize=20)
    

plt.title("2D projection of PSF (RAW Rescaled), E = 5EeV, nside = 128",fontsize=30)
plt.show()
