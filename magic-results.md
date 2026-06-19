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

## 명제 2′ (일반 $M$ — 구조적 분해)

2-진폭 Grover 상태는 항상
$$|\psi\rangle = b\sqrt N\,|{+}\rangle^{\otimes n} + (a-b)\,|\widetilde W\rangle,\qquad
|\widetilde W\rangle=\textstyle\sum_{w\in W}|w\rangle=\sqrt M\,|{\rm flat}_W\rangle,$$
로 쓰이고 $|{+}\rangle^{\otimes n}$ 은 안정자다. 따라서:
- **$W$ 가 아핀부분공간이면** 보조정리 1로 $|{\rm flat}_W\rangle$ 도 안정자 → $|\psi\rangle$ 는
  **안정자 두 개의 중첩** → magic은 **문제크기와 무관하게 유한**(명제 3의 $M{=}1$ 상한 3 bit가
  극단; 수치: 아핀 $M{=}2,4,8$ 정점 $\le3$, $M$ 커지면 감소).
- **$W$ 가 비아핀이면** $|{\rm flat}_W\rangle$ 자체가 magic을 가지며(보조정리 1) 그것이 Grover
  magic에 *추가*된다(수치: 랜덤 $M{=}4$ 정점 $3.29>3$, $M{=}8$ 정점 $4.23$; 같은 $W$의
  ${\rm flat}_W$ magic $=1.54,\,3.36$). **초과 magic의 원천 = 표시집합의 비구조성.**

⟹ 표준 Grover($M{=}1$ 또는 구조적 $W$)의 magic은 유한; 초과 magic은 표시집합 *자체의
비구조성*에서 온다 — 따름정리 1(magic은 (비)구조에 산다)과 같은 주제. *(명시적 닫힌형은
$W$의 자기상관 $A_W(x)=|W\cap(W{\oplus}x)|$에 의존 = 2605.05347의 기하항 $\Lambda$와 동형;
$M{=}1$에서 명제 2로 환원. 일반 닫힌형은 미완.)* `grover_magic.py` §3, `magic_proofs_check.py`.
부호이론적 정정·보정 지표는 `experiments/marker_code_magic.py`: $M_2=-\log_2(\tfrac1{N M^4}\sum_{x,z}\hat g_x(z)^4)$
($g_x(c){=}[c{\in}W\wedge c{\oplus}x{\in}W]$, $\hat{}$는 Walsh), 정정 영점판정 **아핀 $\iff A_W(x){\in}\{0,M\}$**,
비아핀성 스칼라 $\tau$($\tau{=}0\iff M_2{=}0$). 최소 해밍 거리는 지표로 불충분(반례 박제).

---

## 명제 4 (일반 $M$ flat 마커상태 닫힌형 — 가법 에너지)

$|{\rm flat}_W\rangle=|W|^{-1/2}\sum_{x\in W}|x\rangle$ ($M=|W|$, $N=2^n$)에 대해
$$\boxed{\,M_2(|{\rm flat}_W\rangle)=-\log_2\!\Big(\tfrac1{M^4}\sum_{x\in\mathbb F_2^n}E\big(W\cap(W{\oplus}x)\big)\Big)\,},\qquad
E(S)=\#\{(a,b,c,d)\in S^4:a{\oplus}b{\oplus}c{\oplus}d=0\}=\sum_{v}A_S(v)^2,$$
즉 magic은 **이동 자기교집합 $S_x=W\cap(W{\oplus}x)$의 가법 에너지(additive energy)** 합으로 닫힌다.

> **선행 credit(중요):** 이 닫힌형은 **신규가 아니다.** Tarabunga–Castelnovo
> (*Quantum* **8**, 1347 (2024); arXiv:2311.08463) **Eq (8)** 의 임의 상태 4-copy SRE 공식
> $e^{-M_2}=\sum_{\sigma^{(1..4)}}c_{\sigma^{(1)}}c_{\sigma^{(2)}}c_{\sigma^{(3)}}c_{\sigma^{(1)}\sigma^{(2)}\sigma^{(3)}}c^*_{\sigma^{(1)}\sigma^{(2)}\sigma^{(4)}}\cdots$
> 를 균일 flat 상태로 특수화한 것이다(전단사 $(\sigma^{(1)},\sigma^{(2)},\sigma^{(3)},\sigma^{(4)})
> =(a,b,c,a{\oplus}b{\oplus}c{\oplus}x)$ 로 정확히 일치). 본 명제의 가치는 *가법에너지 재명명* +
> 아래 부호이론 특수화(영점판정·Sidon·랜덤 기댓값·Grover 적용)이며, Eq (8) 자체는 2311.08463에 귀속.

*증명(개요).* 실진폭이라 $\langle Z^zX^x\rangle=\tfrac1M\hat g_x(z)$, $g_x=1_{S_x}$.
$\sum_z\hat g_x(z)^4=N\sum_{a\oplus b\oplus c\oplus d=0}1_{S_x}(a)\cdots1_{S_x}(d)=N\,E(S_x)$ (Parseval;
$E(S)=\tfrac1N\sum_z\hat1_S(z)^4$). 대입하면 위 식. **특수값:** $W$ 아핀부분공간이면
$S_x=W\,(x\in V)$ 또는 $\varnothing$, $E(W)=M^3$, $|V|=M$ ⟹ $\sum_x E(S_x)=M^4$ ⟹ $M_2=0$
(보조정리 1 재확인); $M=1$이면 $\sum=E(W)=1$ ⟹ $M_2=0$. 명제 2($M{=}1$)·2′의 미완 닫힌형을
완성. 수치 $\le4\times10^{-15}$ 일치. `experiments/marker_code_closed_form.py`. $\blacksquare$

## 명제 5 (랜덤 마커 통계식 — Sidon 값)

$W$가 **Sidon(B$_2$) 집합**(0 아닌 XOR 차분이 모두 서로 다름)이면 $S_0{=}W$에 $E(W){=}3M^2{-}2M$,
0 아닌 차분 $x$($M(M{-}1)/2$개)마다 $|S_x|{=}2,\ E(S_x){=}8$ 이므로
$$\boxed{\,M_2=\log_2\!\frac{M^3}{7M-6}\ \xrightarrow{M\to\infty}\ 2\log_2 M-\log_2 7\,}.$$
랜덤 $W$는 $M\ll2^{n/2}$에서 whp Sidon이라 $\mathbb E[M_2]\to$ 이 값(수치 $10^{-16}$ 일치). 즉
**구조 없는 마커상태의 magic은 $\sim2\log_2 M$로 자란다** (단 $2\log_2 M$ *증가율*은 표준 한계
$M_\alpha\le2\log_2 R,\ R\le M$ 의 포화로 **기지**; 본 명제의 신규는 *정확 상수 $\log_2 7$·유한형*).
유한 $N$에서 $M^2\gtrsim N$이면 가법
quadruple(충돌)이 늘어 $E\uparrow,\ \xi\uparrow,\ M_2\downarrow$ — 보정은 *하락*이고 $M^2/N$로 통제
(수치: $M{=}8$ gap $0.53\to0.17\to0.00$ as $n{=}6{\to}8{\to}10$). `marker_code_closed_form.py` §2.

## 명제 5′ (랜덤 마커 $\mathbb E[\xi]$ 정확 닫힌형 — $\Delta(M,N)$ 완성)

균일 랜덤 $M$-부분집합 $W\subset\mathbb F_2^n$ ($N{=}2^n$)의 $\xi=2^{-M_2}=\tfrac1{M^4}\sum_x E(S_x)$ 에 대해
$$\boxed{\,\mathbb E[\xi]\,M^4 = \underbrace{(7M^2-6M)}_{\text{Sidon (명제 5)}} + \underbrace{\frac{7\,(M)_4}{N-3} + \frac{N(N-1)(N-2)(N-4)\,(M)_8}{(N)_8}}_{\Delta(M,N)\ (\text{genuine 4중쌍})}\,},\quad (M)_k=\tfrac{M!}{(M-k)!}.$$

*증명(개요).* $\xi M^4=\#\{(a,b,c,d,x):a{\oplus}b{\oplus}c{\oplus}d{=}0,\ \{a,b,c,d\}\cup\{a{\oplus}x,..,d{\oplus}x\}\subseteq W\}$.
base 4중쌍 $N^3$개를 **all-equal**($N$), **paired**($3N(N{-}1)$), **genuine**($N(N{-}1)(N{-}2)$)로 분류
(합 $=N^3$). genuine 4중쌍의 6 차분은 2차원 부분공간 $\{0,\delta_1,\delta_2,\delta_3\}$를 이뤄
$x\in$그 4원소면 $A{\oplus}x{=}A$($k{=}4$), 그 외 $N{-}4$개는 $k{=}8$. 각 8점 배치의 서로 다른 점
수 $k$에 균일 확률 $(M)_k/(N)_k$를 곱해 합산하면 위 식. all-equal+paired가 정확히 Sidon 주항
$7M^2{-}6M$을 준다. **검증:** 닫힌형 vs 몬테카를로 상대오차 $\le1.4\times10^{-2}$(표본오차 수준).
$\blacksquare$

- **점근**($M{\ll}N$): $\Delta\xi\to 7(M)_4/(M^4 N)$ (수치 비$\to1$); $(M)_4/M^4\to1$은 $M\to\infty$.
- $\mathbb E[M_2]\ge-\log_2\mathbb E[\xi]$ (Jensen; $-\log_2$ 볼록), 간극 작아 **$-\log_2\mathbb E[\xi]$ 가
  $\mathbb E[M_2]$의 정밀 추정**. `experiments/marker_code_expected.py`.

## 명제 6 (오라클-가림 = T-비용 가림)

그래프상태 $|\psi_f\rangle=2^{-n/2}\sum_x|x\rangle|f(x)\rangle$ 와 오라클 $U_f:|x\rangle|y\rangle\mapsto|x\rangle|y{\oplus}f(x)\rangle$ 에 대해
$$\boxed{\,M_2(|\psi_f\rangle)>0\iff f\text{ 의 ANF에 차수}\ge2\text{ 단항식 존재}\iff U_f\text{ 게이트분해가 비클리포드(Toffoli/}T)\text{ 게이트 요구}\,}.$$
세 조건 모두 **$f$ 아핀에서만 동시 0**(따름정리 1 + ANF: 1차항=CNOT 클리포드, 2차항=Toffoli=$T$).

*근거.* 따름정리 1로 $M_2{=}0\iff f$ 아핀$\iff$ ANF가 1차 이하. 출력비트 ANF의 차수$\ge2$ 단항식
수 $T_{\rm proxy}$(=Toffoli/$T$ 비용 대용)와 $M_2$가 **함께 0에서 켜지고 함께 증가**(수치 $n{=}4$:
$T_{\rm proxy}{=}0,1,2,3,4 \Rightarrow M_2{=}0,1.54,2.48,3.70,4.43$). 끝점: **Simon**(선형 오라클,
$T_{\rm proxy}{=}0,\ M_2{=}0$인데도 지수 쿼리 우위) vs **Shor** modexp(비선형, $T_{\rm proxy}{=}4{\to}156,\
M_2{=}1.54{\to}5.23$). ⟹ 쿼리 오라클이 숨기는 "보이지 않는 magic"은 곧 게이트분해의 $T$-비용이다.
*(2507.16543의 Clifford-가림과 구분: 이쪽은 오라클/추상화가 가림.)* `experiments/oracle_tcount_magic.py`.
**FTQC 환산(`experiments/oracle_ftqc_estimate.py`):** 차수 $d{\ge}2$ ANF 단항식 $\to(2d{-}3)$ Toffoli
$\to7T$(ancilla $4T$). $T_{\rm est}{=}0\iff M_2{=}0$; modexp ANF-합성 상한 $N{=}15{\to}4$ Toffoli
$\sim N{=}35{\to}4171$. **단 출력비트별 ANF 합성은 상한**(산술 오라클은 windowed로 훨씬 낮음 — 실
Shor는 Gidney–Ekerå 2021). 실용: 정확한 영점판정 + 정성 magic↔T 법칙 + 구조없는 오라클 상한 지표.

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
