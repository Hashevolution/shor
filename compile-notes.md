# LaTeX 컴파일·arXiv 업로드 가이드

## 로컬 컴파일

LaTeX 배포판 (Windows 권장: MiKTeX, macOS: MacTeX, Linux: TeX Live) 설치 후:

```powershell
pdflatex paper.tex
pdflatex paper.tex   # 두 번 실행 — 참조 해석 (cross-references, citations)
```

또는 `latexmk` 사용:

```powershell
latexmk -pdf paper.tex
```

산출물: `paper.pdf`.

## 컴파일 시 자주 빠지는 패키지

`paper.tex` 가 의존하는 패키지 (MiKTeX 은 자동 다운로드, TeX Live 는 사전 설치 필요):

- `amsmath`, `amssymb`, `amsthm` — 수식·정리
- `geometry` — 마진
- `hyperref` — 링크
- `booktabs`, `array` — 표
- `listings`, `xcolor` — 코드 블록
- `microtype` — 조판 미세 조정
- `babel` (english) — 영어

## TODO before arXiv 업로드

1. **저자 정보 채우기**. `paper.tex` 의 `\author{Author Name}` 와 `\thanks{TODO: ...}` 를 실명·소속·이메일로 교체.
2. **날짜 확인**. `\date{June 2026}` — 필요시 갱신.
3. **인용 누락 확인**. `pdflatex` 출력에서 `Warning: Citation ... undefined` 가 없는지.
4. **그림 없음 확인**. 본 paper 는 표만 사용 — 그림 누락 경고 없음.
5. **arXiv 카테고리 결정**. 권장: `quant-ph` (Quantum Physics) 주 + `cs.CC` (Computational Complexity) 또는 `math.NT` (Number Theory) 보조.
6. **§3.6 의 최신 결과 통합 확인**. v0.2 (13-seed × 12 σ + (1147, 2) cross-cell + amplification) 모두 반영되었는지.

## Zenodo 업로드 (arXiv 와 *별도* 가능)

arXiv endorsement 대기 중에 *Zenodo 만* 으로 priority lock 가능:

1. https://zenodo.org/ 계정 + GitHub OAuth 연결
2. https://zenodo.org/account/settings/github/ 에서 `Hashevolution/shor` 토글 ON
3. GitHub repo PUBLIC 으로 전환
4. GitHub Release 생성:
   ```powershell
   gh release create v0.2.0 --title "Multi-boundary mechanism observation" `
       --notes-file release_notes_v0.2.0.md
   ```
5. Zenodo 가 자동 archive + DOI 발급 (~10분)
6. CITATION.cff / README.md 에 DOI badge 추가

## arXiv 업로드 단계

1. https://arxiv.org/submit 접속 (계정 필요, endorsement 필요할 수 있음 — 첫 quant-ph 게재면 누군가의 보증 필요)
2. "Start submission" → "Article" 선택
3. **소스 업로드**: `paper.tex` 단독 (그림·BibTeX 없음 → 단일 파일 충분).
   - arXiv 가 자동으로 컴파일 수행. 실패 시 로그 확인.
4. **메타데이터 입력**:
   - Title: "A Noise-Invariant Determinism Theorem for Multi-Base Post-Processing in Shor's Order Finding"
   - Authors: (저자 정보)
   - Abstract: paper.tex 의 abstract 복사
   - Categories: `quant-ph` (primary), `cs.CC` 또는 `math.NT` (secondary)
   - MSC class (선택): `81P68` (Quantum computation) 또는 `11Y05` (Factorization)
   - ACM class (선택): `F.2.1` (Numerical Algorithms and Problems)
   - Comments: "8 pages. Companion code: https://github.com/Hashevolution/shor"
5. **License**: 권장 CC BY 4.0 또는 arXiv perpetual non-exclusive.
6. **Preview** → **Submit**. 게재 승인 대기 (보통 1 영업일).

## arXiv 거절 케이스 (사전 자체 점검)

- 너무 짧음 (< 4 페이지): paper 약 8 페이지 → ok
- 영어 품질: 수동 한 번 읽어보기 권장
- 기존 발표·게재물 중복: 없음
- 학술 콘텐츠 부족: 본 paper 는 정리 + 검증 + 음수 결과 + 인용 — 충분
- 익명 저자: 본명 필수

## 엠바고 / 후속

- 게재 후 paper 가 검색 가능 (arXiv ID 부여)
- 인용 시: `arXiv:YYMM.NNNNN`
- 후속 버전: arXiv 가 v2, v3 ... 자동 관리

## 회의/저널 후속 투고

본 paper 는 다음 venue 에 적합 수준:

- **Workshop**: TQC (Theory of Quantum Computation), QIP (workshop track), AQIS
- **저널**: Quantum Information Processing (Springer), Quantum Information \& Computation (Rinton)
- **단편 (short note)**: arXiv 단독 게재가 가장 현실적
