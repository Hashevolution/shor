"""
magic.py — 순수상태의 stabilizer 2-Rényi 엔트로피(비안정자성, "magic") 측정.

정의 (Leone–Oliviero–Hamma, PRL 128 050402):
    M₂(|ψ⟩) = -log₂( (1/2ⁿ) Σ_{P∈𝒫ₙ} ⟨ψ|P|ψ⟩⁴ ),   안정자 상태 ⟹ M₂ = 0,  0 ≤ M₂ ≤ n.

전수합은 4ⁿ개 Pauli지만, 다음 항등식으로 XOR-Walsh–Hadamard 변환을 써서 O(n·4ⁿ)에
정확 계산한다 (cf. arXiv:2512.24685 XOR-FWHT; arXiv:2601.07824 state-vector magic):

    ⟨ψ|ZᶻXˣ|ψ⟩ = (-1)^{z·x} · WHT_z[ h_x ](z),   h_x(c) = ψ*_{c⊕x} ψ_c,
    ⟨P⟩ 는 실수(Hermitian Pauli)이므로  Σ_P ⟨P⟩⁴ = Σ_{x,z} |WHT[h_x](z)|⁴.

`sre2` 가 주 도구이고, `sre2_bruteforce` 는 Pauli 행렬을 직접 만들어 검산한다(n ≤ 5).
이 모듈은 Grover/Simon/Shor comb 상태 어디에도 그대로 쓰인다.
"""
from __future__ import annotations

import itertools

import numpy as np


def _fwht_1d(a: np.ndarray) -> np.ndarray:
    """Unnormalized Walsh–Hadamard transform of a length-2ⁿ vector (complex ok)."""
    a = np.asarray(a, dtype=complex).copy()
    N = a.size
    h = 1
    while h < N:
        a = a.reshape(N // (2 * h), 2, h)
        x = a[:, 0, :]
        y = a[:, 1, :]
        a = np.stack([x + y, x - y], axis=1).reshape(N)
        h *= 2
    return a


def sre2(psi: np.ndarray) -> float:
    """순수상태 |ψ⟩(길이 2ⁿ 진폭벡터)의 stabilizer 2-Rényi 엔트로피 M₂ (bits)."""
    psi = np.asarray(psi, dtype=complex).ravel()
    N = psi.size
    n = int(round(np.log2(N)))
    if 2 ** n != N:
        raise ValueError("state length must be a power of two")
    nrm = np.linalg.norm(psi)
    if nrm == 0:
        raise ValueError("zero state")
    psi = psi / nrm

    idx = np.arange(N)
    total = 0.0
    for x in range(N):
        hx = np.conj(psi[idx ^ x]) * psi          # h_x(c) = ψ*_{c⊕x} ψ_c
        W = _fwht_1d(hx)
        total += np.sum(np.abs(W) ** 4)           # Σ_z |WHT[h_x](z)|⁴
    xi = total / N                                 # (1/2ⁿ) Σ_P ⟨P⟩⁴
    m2 = float(-np.log2(xi))
    return 0.0 if abs(m2) < 1e-9 else m2            # M₂≥0; FP 잡음(±0 포함)만 0으로 정리


# ── 검산용 (느림, n ≤ 5) ─────────────────────────────────────────────────────
_P1 = {
    0: np.array([[1, 0], [0, 1]], dtype=complex),       # I
    1: np.array([[0, 1], [1, 0]], dtype=complex),       # X
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),    # Y
    3: np.array([[1, 0], [0, -1]], dtype=complex),      # Z
}


def sre2_bruteforce(psi: np.ndarray) -> float:
    """4ⁿ Pauli를 직접 만들어 M₂를 계산 (sre2 검산용)."""
    psi = np.asarray(psi, dtype=complex).ravel()
    N = psi.size
    n = int(round(np.log2(N)))
    psi = psi / np.linalg.norm(psi)
    total = 0.0
    for combo in itertools.product(range(4), repeat=n):
        P = np.array([[1.0 + 0j]])
        for p in combo:
            P = np.kron(P, _P1[p])
        ev = (psi.conj() @ (P @ psi)).real
        total += ev ** 4
    return float(-np.log2(total / N))
