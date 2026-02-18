# Analysis Scripts

`/home/hea/Kyle/Analysis`에서 분석 스크립트만 추려 정리한 저장소입니다.

## Included
- `scripts/`: 메인 분석 스크립트
- `scripts/testscript/`: 보조 테스트/시각화 스크립트

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
