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

- `CutEnergy.py`: This returns the Cut5EeV.csv file from Archieve_v6r2p2 file.
- D

## Excluded
- 대용량 데이터(`*.csv`)
- 결과 이미지(`*.png`)
- 백업 파일(`*~`, `#*#`)
- 아카이브 파일(`Archive_*`)

## Setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Notes
- 스크립트는 원본 파일명을 유지했습니다.
- 실행에 필요한 입력 데이터(`.csv`)는 별도 경로에 두고 파일명/경로를 스크립트에 맞게 설정해야 합니다.
