# Magic과 속도우위 — 명제와 증명 (T3)

*`magic-and-quantum-speedup.md`의 수치 결과를 명제로 정리. 측정도구 `magic.py`,
재현·검증 `experiments/magic_proofs_check.py`(모든 명제를 작은 $n$에서 assert).*

표기: $n$큐비트 순수상태 $|\psi\rangle$의 stabilizer 2-Rényi 엔트로피
$$M_2(|\psi\rangle)=-\log_2\!\Big(\tfrac1{2^n}\sum_{P\in\mathcal P_n}\langle\psi|P|\psi\rangle^4\Big),\qquad 0\le M_2\le n,\quad M_2=0\iff |\psi\rangle\ \text{안정자}.$$

---

## 보조정리 1 (평탄상태의 안정자성)

받침 $S\subseteq\mathbb F_2^n$ 위에 균등진폭·동일위상인 상태
$|\psi_S\rangle=|S|^{-1/2}\sum_{x\in S}|x\rangle$ 에 대해
$$\boxed{\,M_2(|\psi_S\rangle)=0 \iff S\ \text{가}\ \mathbb F_2^n\ \text{의 아핀부분공간}\,}$$
(즉 어떤 선형부분공간 $V$와 $v$에 대해 $S=v\oplus V$).

*증명.* ($\Leftarrow$) $S=v\oplus V$이면 $|\psi_S\rangle$은 $V$의 $X$형 생성원
$\{X^{u}:u\in V\}$의 적절한 부분과 $V^\perp$의 $Z$형 체크로 안정화되는 안정자상태이므로
$M_2=0$. ($\Rightarrow$) 안정자 순수상태의 받침은 항상 아핀부분공간이며, 그 위 진폭은
$2^{-k/2}i^{\ell(x)}(-1)^{q(x)}$ 꼴($\ell$ 선형, $q$ 2차)이다 (Dehaene–De Moor 2003;
Gross 2006). $M_2=0\iff$ 안정자상태이고, 균등진폭·동일위상은 $\ell=q=0$인 경우이므로
받침이 아핀부분공간이다. $\qquad\blacksquare$

---

## 따름정리 1 (오라클: Simon vs Shor)

함수 $f:\mathbb F_2^n\to\mathbb F_2^m$의 그래프상태
$|\psi_f\rangle=2^{-n/2}\sum_x|x\rangle|f(x)\rangle$ 에 대해
$$\boxed{\,M_2(|\psi_f\rangle)=0 \iff f\ \text{가}\ \mathbb F_2\text{-아핀}\ (f(x)=Ax\oplus b)\,}.$$

*증명.* 받침은 $\mathrm{graph}(f)=\{(x,f(x)):x\}\subseteq\mathbb F_2^{n+m}$. 이것이
아핀부분공간 $\iff$ $(x,f(x))$가 $x$의 아핀함수 $\iff$ $f$ 아핀. 보조정리 1. $\blacksquare$

- **Simon.** 숨은문자열 $s$의 2-대-1 함수를 선형으로 잡을 수 있다:
  $f(x)=x\oplus(x_p\!\cdot\! s)$ ($p$는 $s$의 set bit), $\ker=\{0,s\}$, 선형. 따라서
  오라클이 클리포드 → **회로 전체 $M_2\equiv0$이면서 지수 쿼리 속도우위.** magic은 (쿼리모델)
  지수 속도우위에 *필요하지 않다.*
- **Shor.** $f(x)=a^x\bmod N$ 은 곱셈적=비선형 → $M_2>0$ (수치: $N{=}15{\to}1.54$,
  $21{\to}5.80$, $33{\to}7.12$). magic은 modexp의 **비선형성**에서 나온다.
- **함의(오라클-가림).** magic의 원천은 회로 표면이 아니라 *문제의 비선형성*이며, 쿼리 오라클
  (과 `shor.py`의 `np.fft` 지름길)이 그 magic을 블랙박스로 *숨긴다*.

## 따름정리 2 (측정후 comb)

$\mathrm{comb}_r=\{x\in\mathbb Z_{2^t}: x\equiv x_0\ (\mathrm{mod}\ r)\}$ 위 균등중첩에 대해
$M_2=0\iff$ 그 잉여류가 $\mathbb F_2^t$의 아핀부분공간.
- $r=2^s$ ⟹ 하위 $s$비트 고정 = 아핀부분공간 ⟹ $M_2=0$.
- $r$에 홀수 인수 ⟹ 일반적으로 아핀 아님 ⟹ $M_2>0$.
- *예외(우연한 선형화):* $r(m{-}1)=2^t{-}1$ 류는 받침이 $\{(j{\ll}s)\,|\,j\}$ 꼴 선형부분공간이 되어
  $M_2=0$ (예 $t{=}10,r{=}33$: $33{\cdot}31{=}1023$). 따라서 "$M_2{=}0\iff r$이 2의 거듭제곱"은
  *거의*-참이고 정확한 판정은 보조정리 1(받침의 아핀성).

---

## 명제 2 (Grover 단일표시 닫힌형)

표시 1개에 진폭 $a$, 나머지 $N{-}1$개에 $b$, $a^2+(N{-}1)b^2=1$인 상태($N=2^n$)에 대해
$$\sum_{P}\langle P\rangle^4 = 1+(N{-}1)(a^2{-}b^2)^4+(N{-}1)\big((N{-}2)b^2{+}2ab\big)^4+(N{-}1)\big(\tfrac N2{-}1\big)\big(2b(a{-}b)\big)^4,$$
$$M_2=-\log_2\!\Big(\tfrac1N\sum_P\langle P\rangle^4\Big).$$

*증명(개요).* 실진폭이라 $\langle Z^zX^x\rangle=\sum_c\psi_c\psi_{c\oplus x}(-1)^{z\cdot c}$, 그리고
$i$위상은 4제곱에서 사라져 $\sum_P\langle P\rangle^4=\sum_{x,z}\langle Z^zX^x\rangle^4$. 표시를 0번지로
두고 $x{=}0$/$x{\ne}0$, $z{=}0$/$z{\cdot}x{=}0$/$z{\cdot}x{=}1$로 분류해 직접 합산하면 위 식.
($z{\cdot}x{=}1,x{\ne}0$ 항은 0.) 수치와 $10^{-13}$ 일치. $\blacksquare$

## 명제 3 (Grover 정점 magic $\to 3$ bit)

명제 2에서 $N\to\infty$, $a$ 고정($b^2\!\approx\!(1{-}a^2)/N$)이면
$$M_2 \to -\log_2\!\big(a^8+(1-a^2)^4\big),$$
이는 $a^2=\tfrac12$에서 최대 $-\log_2(1/8)=\boxed{3\ \text{bit}}$. 따라서
$$\sup_k M_2(\text{Grover}_n) \xrightarrow[n\to\infty]{} 3,\qquad \text{밀도}\ M_2/n\to 0.$$

*증명.* $b^2(N{-}2)\to1{-}a^2$, $2ab\to0$이므로 명제 2의 $S/N\to a^8+(1-a^2)^4+O(1/N)$.
$g(u)=u^4+(1-u)^4$ ($u=a^2$)는 $g'(u)=4u^3-4(1-u)^3=0\Rightarrow u=\tfrac12$에서 최소
$g=\tfrac18$. 그러므로 $M_2$의 상한 $\to-\log_2\tfrac18=3$. $\blacksquare$

*경계.* 명제 2–3은 단일표시($M{=}1$)에 대한 것. 일반 $M$은 표시집합 $W$의 자기상관
$A_W(x)=|W\cap(W{\oplus}x)|$에 의존(2605.05347의 기하항 $\Lambda$와 동형)하며, 수치상
고정 $M$의 정점은 여전히 $O(1)$(밀도 낮음, `grover_magic.py` §3). 해석적 일반화는 미완.

---

## 속도우위 사다리 (요약)

| 알고리즘 | 속도우위 | $M_2$ 거동 | 밀도 $M_2/n$ | 근거 |
|---|---|---|---|---|
| Simon / BV | 지수(쿼리) | $0$ (아핀 오라클=클리포드) | $0$ | 따름정리 1 |
| Grover | 2차(다항) | 정점 $\to 3$ bit, 답서 $0$ | $\to 0$ | 명제 2–3 |
| Shor (comb) | 지수 | $\propto t$ 증가 | $\sim0.4$–$0.55$ | 따름정리 2 |
| Shor (in-circuit) | 지수 | $\to L$ (최대) | $\to 1$ | 2605.05347 |

> **결론.** 비안정자성의 *양/밀도*가 속도우위의 유형을 가른다. magic은 회로 표면이 아니라
> 문제의 비선형성(아핀 vs 곱셈적)에 살며, 오라클·FFT 같은 블랙박스가 그것을 숨길 수 있다.
> Grover의 2차 속도우위는 *유한한* magic을 썼다가 되돌리는 반면, Shor의 지수 속도우위는
> magic을 문제크기와 함께 키운다.

*인접 선행과의 선:* 양자걷기 magic(2506.17783·2504.19750)은 1D 격자 수송으로 *시간축*
포화이며 완전그래프/탐색을 다루지 않음(전문 확인). 알고리즘 magic 비교(2505.17185·2507.16543)는
변분/QFT 회로이고 Grover·Simon·Shor 사다리는 다루지 않음. 위 명제 2–3·따름정리 1은 그 빈칸.
