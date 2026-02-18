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

### 
