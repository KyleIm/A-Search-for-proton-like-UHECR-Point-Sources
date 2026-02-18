import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from crpropa import *
import pandas as pd

df=pd.read_csv('Cut50EeV.csv')

l=df['L'].values * -1
b=df['B'].values


#print(dL)

plt.figure()

ax = plt.subplot(111, projection='mollweide')
ax.grid(True)
ax.scatter(l, b, s=10)
xticks = np.radians([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
ax.set_xticks(xticks)
ax.set_xticklabels(["180°", "150°", "120°", "90°", "60°", "30°", "0°", "-30°", "-60°", "-90°", "-120°", "-150°", "-180°"], fontsize=20)
#xticks = ax.get_xticks()
#ax.set_xticklabels([f"{-int(np.degrees(tick))}°" for tick in xticks], fontsize=20)

ax.tick_params(axis='y', labelsize=20)  # Y-axis tick font size

# Add axis labels and legend
ax.set_xlabel("L", fontsize=25)
ax.set_ylabel("B", fontsize=25)
#ax.legend(fontsize=25)
plt.title("Selected Event location for magnetic field", fontsize=30)
plt.show()
