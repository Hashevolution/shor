"""
대형 N 에서 (C)-determinism 정리 검증 실행 + 결과 파일 저장.
재실행 시 `python verify_large_run.py` 로 paper §4 확장표 재현.
"""

import sys
import time
from contextlib import redirect_stdout

from demo import verify_c_determinism


SUBSET = ["noise-free", "depol p=0.8", "phase σ=2.5", "modexp q=0.8"]
TRIALS = 100
DEFAULT_NS = [1147, 2491, 4087]


def main(argv):
    Ns = [int(x) for x in argv[1:]] if len(argv) > 1 else DEFAULT_NS
    out_path = "verify_large_results.txt"
    lines = []

    class Tee:
        def write(self, s):
            sys.__stdout__.write(s)
            lines.append(s)
        def flush(self):
            sys.__stdout__.flush()

    tee = Tee()
    with redirect_stdout(tee):
        print(f"=== verify_large: N={Ns}, trials={TRIALS}, noise={SUBSET} ===")
        for N in Ns:
            t0 = time.time()
            verify_c_determinism(N, trials=TRIALS, noise_subset=SUBSET)
            print(f"  [N={N} elapsed: {time.time()-t0:.1f}s]\n")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("".join(lines))
    print(f"\n저장: {out_path}", file=sys.__stdout__)


if __name__ == "__main__":
    main(sys.argv)
