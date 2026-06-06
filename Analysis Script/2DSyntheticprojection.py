import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp
from healpy.newvisufunc import projview

df=pd.read_csv('Sy_events_5EeV.csv')
l=df['L'].values
b=np.pi/2-df['B'].values
dl=df['dL'].values
db=np.pi/2-df['dB'].values


nside=15
npix=hp.nside2npix(nside)

indices=hp.ang2pix(nside,b,l)
d_indices=hp.ang2pix(nside,db,dl)

hpmap=np.zeros(npix)

for i in range(len(l)):
    hpmap[indices[i]]=hpmap[indices[i]]+1

hpmap_masked = np.ma.masked_where(hpmap == 0, hpmap)

projview(
    hpmap_masked,
    coord=["G"],
    graticule=True,
    graticule_labels=True,
    xlabel="l",
    ylabel="b",
    flip="astro",
    projection_type="mollweide",
)

ax = plt.gca()

# Mark selected HEALPix pixel centers with red "X" markers.
target_indices = np.array([2100, 2150, 2200], dtype=int)
theta_c, phi_c = hp.pix2ang(nside, target_indices)
lat_c = np.pi / 2.0 - theta_c
lon_c = (phi_c + np.pi) % (2.0 * np.pi) - np.pi
lon_c = -lon_c  # flip="astro"
marker_colors = ["red", "deepskyblue", "lime"]
marker_size = 176  # 20% smaller than 220
counts = np.bincount(d_indices, minlength=npix)

for i, pix in enumerate(target_indices):
    ax.scatter(
        lon_c[i],
        lat_c[i],
        marker="x",
        s=marker_size,
        c=marker_colors[i],
        linewidths=3.0,
        zorder=10,
        label=f"Indices {pix} : {int(counts[pix])}",
    )

print("Counts from dL,dB-based HEALPix indices:")
for pix in target_indices:
    print(f"{pix}: {int(counts[pix])}")

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(20)  

ax.set_xlabel("L", fontsize=25)
ax.set_ylabel("B", fontsize=25)
    
cbar = plt.gcf().axes[-1]
cbar.tick_params(labelsize=20)
plt.title(r"2D projection of Synthetic Events, E > 5EeV, $N_{side}$ = 15",fontsize=30)

legend_handles, legend_labels = ax.get_legend_handles_labels()
ax.legend(
    legend_handles,
    legend_labels,
    loc="upper right",
    fontsize=20,
    framealpha=0.9,
)


plt.show()
