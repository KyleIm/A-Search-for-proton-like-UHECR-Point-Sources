# Kyle Im's PhD analysis Repository

Welcome to the repository, which includes all the important source codes that Kyle used during his PhD analysis.

## Introduction
This repository includes 3 main parts.

### Analysis Script 
Here you can find all the analysis scripts that Kyle has used. Users can reproduce Kyle's result just by running the same Python script saved here. However, due to the storage limitations of this repository, users are advised to obtain the Archive_v6r2p2 files directly from the Auger Herald, as they constitute the core dataset for this analysis. Please type 
```
wget http://physics-anduril.case.edu/~covault/herald/Archive_v6r2p2
```
for download. The file will be approximately 1.6 Gb.

Also, for the continuity for the setup, when we use the Energy, we cut the energy for the 5 EeV threshold. You can change the energy preference based on what you want, however, we will explain based on this setup.

Followings are the list of scripts saved in this folder and what they do.

- `CutEnergy.py`: This returns the Cut5EeV.csv file from Archieve_v6r2p2 file. The Cut5EeV.csv file is a smaller capacity file that only includes the necessary information for the analysis.
- `Deflection.py`: This returns the JF12_5EeVCut.csv file from Cut5EeV.csv file. JF12_5EeVCut.csv includes JF12 Backtracked coordinate information. Additional parameters may be introduced depending on the direction of the simulation. (However, we did not introduce it in this analysis.)
- `TimeShuffling.py`: This returns test5EeV.csv file from Cut5EeV.csv file. The test5EeV.csv file includes simulated dataset from time shuffling simulation.
- `location_selection.py` : This creates a scattering projection plot using Mollewide projection. We are using the astro convention, so the galactic longitude is increasing from right to left. Since the default setup for Geographic Projections in matplotlib is increasing from left to right, some modifications were made for this.
- `RPS_5EeV.py` : This creates a scattering on the Mollewide projection plot for random parameter scattering. The default setup is a backtracking from the origin with a 5 EeV energy test particle with 1,000 samplings. Randomized JF12 sub-parameters make different results whenever we make the backtracking, and this is represented as blue. For comparison, we have represented a green dot for the origin of the particle and a red dot for the backtracking result using JF12 model. The result is shown as  `RPS.png`, however, the plot on the thesis is a modified version that uses the name `Deflection_from_origin.png`. (Modification made with google slide.)
- `Firstordercandidate.py` : This shows the projection for selected locations of 221 Pierre Auger events above 50 EeV. These events are selected for the purpose of choosing the location where we create random parameter sampling.
- `combined_evrs_effectivearea.py` : This script performs 3 objectives. First, it shows Energy vs Average Distance plot. Here Average Distance stands for average distance between true value (when JF12 does not have uncertainty) and simulated value (when we apply uncertainty for JF12) for 221,000 runs. (1,000 runs for 221 selected locations) and add dotted guideline what would be the equivalent HealPix N_{side} shows the distribution of effective distances at TARGER_ENERGY, shown on a log scale (y-axis). We use 221 locations which we have chosen from `Firstordercandidate.py`. The default value for TARGER_ENERGY would be 5 EeV.

### Modification for CRPropa3

Customization for CRPropa3 was made for the purpose of Gaussian random parameter sampling for JF12 sub parameters. When you download the CRPropa3, please replace `JF12Field.h` and `JF12Field.cpp` file from 
```
CRPropa3/src/magneticfield/
```
before compiling the library. After implementing customization, you can use Gaussian random sampling for JF12 magnetic field model when you use JF12Field(True). The default JF12Field() is False, so if you do not want the random sampling, you can still use JF12Field() as written on CRPropa3 manual.

### Updated plots from Kyle's Thesis

This folder contians most updated plots from the latest version of Kyle's analysis plot. The main purpose for this folder is to update  includes all the plot which Kyle made for the Chapter 3 and 4
