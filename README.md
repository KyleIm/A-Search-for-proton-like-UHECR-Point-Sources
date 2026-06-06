# Kyle Im's PhD analysis Repository

Welcome to the repository, which includes all the important source codes that Kyle used during his PhD analysis.

## Introduction
This repository includes 4 main parts. Analysis Scripts, CSV Files, Modifications for CRPropa3, Updated Plots for Kyle's Thesis.

### Analysis Scripts 
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
- `location_selection.py` : This shows the projection for selected locations of 221 Pierre Auger events above 50 EeV using Mollewide projection. We are using the astro convention, so the galactic longitude is increasing from right to left. Since the default setup for Geographic Projections in matplotlib is increasing from left to right, some modifications were made for this.
- `RPS_5EeV.py` : This creates a scattering on the Mollweide projection plot for random parameter scattering. The default setup is a backtracking from the origin with a 5 EeV energy test particle with 1,000 samplings. Randomized JF12 sub-parameters make different results whenever we make the backtracking, and this is represented as blue. For comparison, we have represented a green dot for the origin of the particle and a red dot for the backtracking result using JF12 model. The result is shown as  `RPS.png`, however, the plot on the thesis is a modified version that uses the name `Deflection_from_origin.png`. (Modification made with google slide.)
- `combined_evrs_effectivearea.py` : This script performs 3 objectives. First, it shows histograms for TARGET_ENERGY (the default setup is 5 EeV) for average distance. Here, the Average Distance stands for the average distance between the true value (when JF12 does not have uncertainty) and simulated value (when we apply uncertainty for JF12) for 221,000 runs. (1,000 runs for 221 selected locations) It shows plot in two formats. One is semilog, and the other is log-log. Second, it shows Energy vs Average Distance plot for different energy (E from 1 EeV to 10 EeV. Increases from 0.5 EeV). Also it addes dashed guideline what would be the equivalent area for given HealPix N<sub> side </sub>. Lastly, it saves information of E, mean, sigma, nearest_nside, nearest_nside_value as `EvsR.txt` and `EvsR.csv` file.
- `2Dprojection.py` : This shows you the HealPix Mollweide projection for the cosmic ray event which you want to project. This is basic version of projection plot, so there are variation of this script such as `lima_healpyprojection.py`, `Flux_projection.py`, `2DSyntheticprojection.py`, or lima_healpyprojection_Sy.py. You can project RAW dataset (Use the csv file which you obtained from `CutEnergy.py`) or JF12 backtracked dataset (Use the csv file which you obtained from `Deflection.py`) or use the simulated dataset from time shuffling. (Use the csv file which you obtained from `TimeShuffling.py`) Default is nside = 15, E = 5 EeV, with file named `JF12_5EeVCut.csv`. However, you can change this setup based on what you want to project.
-  `lima_healpyprojection.py` : This is a variation of `2Dprojection.py`. This shows HealPix projection of Li-Ma significance. Default setting is for E = 5 EeV and Li-Ma significance for JF12 Backtracked events. You can change this with a different filename and by changing the setting. However, be sure to rescale the expected number of events or simulations to be the same.
-  `lima_histogram.py` : This shows histogram of Li-Ma significance for both RAW and JF12 backtracked events. Dashed Gaussian fitting lines and the textbox for the result of the plots are also applied to the histogram. When you run this script, it could show an error warning when there is a bin which have real data but without simulated data. In this case, the bin is excluded.
-  `Flux_projection.py` : This is a variation of `2Dprojection.py`. This shows a HealPix projection of the flux upper limit. The hardware acceptance gets maximized at E > 3 EeV. Not able to use this if the energy threshold is lower than this. The default energy setting is for E > 5 EeV. If you are using a different energy scale projection, I recommend renormalize scale based on the number of events and number of HealPix bins. (Rescale the factor of expectation of number of event or simulations to be the same.)
-   `Flux_projection_low.py` : This is a variation of `Flux_projection.py`. The energy threshold for this script is 2 EeV. We have directional exposure for this lower energy threshold.
-   `Decvsflux.py` : This shows declination dependence of the flux upper limit as a step plot for events with a 5 EeV threshold. I have made an average flux for each step, and the step is made with a 3-degree order.
-   `Decvsflux_all.py` : This is a variation of `Decvsflux.py`. Now this shows step functions for 3 different energy thresholds. E > 5 EeV is shown as red, E > 3 EeV is shown as blue and 2 EeV < E $\leq$ 3 EeV is shown as green. The flux is rescaled based on 5 EeV.
-   `lima_Syenthetic.py` : This creates `Sy_events_5EeV.csv` which used for `2DSyntheticprojection.py`, `lima_healpyprojection_Sy.py`, and `lima_histogram_Syn.py`. This csv file is synthetic events used for synthetic event injection test.
-   `2DSyntheticprojection.py` : This is a variation of `2Dprojection.py`. The purpose for this projection plot is to indicate 3 bins where we apply synthetic events. This requires `Sy_events_5EeV.csv` file.
-   `lima_healpyprojection_Sy.py` : This is a variation of `2Dprojection.py`. The projection is produced from the sum of the outputs generated by `2Dprojection.py` and `2DSyntheticprojection.py`. This requires `Sy_events_5EeV.csv` file.
-   `lima_histogram_Syn.py` : This is a variation of `lima_histogram.py`. This shows a histogram including synthetic events. This requires `Sy_events_5EeV.csv` file.

### CSV Files

The purpose of 

### Modifications for CRPropa3

Customization for CRPropa3 was made for the purpose of Gaussian random parameter sampling for JF12 sub parameters. When you download the CRPropa3, please replace `JF12Field.h` and `JF12Field.cpp` file from 
```
CRPropa3/src/magneticfield/
```
before compiling the library. After implementing customization, you can use Gaussian random sampling for JF12 magnetic field model when you use JF12Field(True). The default JF12Field() is False, so if you do not want the random sampling, you can still use JF12Field() as written on CRPropa3 manual.

### Updated plots from Kyle's Thesis

This folder contians most updated plots from the latest version of Kyle's analysis plot. The main purpose for this folder is to update  includes all the plot which Kyle made for the Chapter 3 and 4
