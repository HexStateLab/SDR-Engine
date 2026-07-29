#!/usr/bin/env python3
"""
room_find_a.py — Deterministic equation for the room's Shor's algorithm.

Given N (integer to factor) and Q (register size, default=256×banks),
computes the smallest coprime a whose period r = period(a,N) satisfies:
  1. r is even
  2. Some QFT peak k*Q/r produces a continued-fraction convergent
     denominator q where q × m = r for m ∈ {1,2,3,4,8,16}.

This a is the value the room's EM-field Shor's will successfully factor N with.

Usage:
  python3 room_find_a.py <N> [banks] [--all]
    N      — integer to factor
    banks  — number of frequency pages (default 1, Q=256)
    --all  — show all compatible (a,r) pairs, not just the first

Examples:
  python3 room_find_a.py 2021 32      # a for N=2021 with 32 banks
  python3 room_find_a.py 3127 64      # a for N=3127 with 64 banks
  python3 room_find_a.py 77           # a for N=77 with default Q=256
"""

import math
import sys


def period(a, N):
    """Find the multiplicative order of a modulo N (period of a^x mod N)."""
    seen = {}
    x = 1
    for i in range(N + 2):
        seen[x] = i
        x = (x * a) % N
        if x in seen:
            return i + 1 - seen[x]
    return -1


def continued_fraction_denominators(num, den):
    """Return all denominator values from the continued fraction of num/den."""
    a, b = num, den
    coeffs = []
    while b:
        coeffs.append(a // b)
        a, b = b, a % b

    if not coeffs:
        return []

    denominators = []
    p0, q0 = 0, 1
    p1, q1 = 1, coeffs[0]
    denominators.append(q1)

    for i in range(1, len(coeffs)):
        p = coeffs[i] * p1 + p0
        q = coeffs[i] * q1 + q0
        denominators.append(q)
        p0, q0 = p1, q1
        p1, q1 = p, q

    return denominators


def room_compatible(N, a, Q=256):
    """
    Check if (N,a) is compatible with the room's Shor's at register size Q.
    Returns (True, r, c, q, m) or (False, r, None).
      r = period(a,N)
      c = QFT peak bin position
      q = continued-fraction denominator
      m = multiplier such that q×m = r
    """
    r = period(a, N)
    if r < 2 or r % 2 == 1:
        return False, r, None

    for k in range(1, min(r, 500)):
        c = round(k * Q / r)
        if c <= 0 or c >= Q:
            continue
        g = math.gcd(c, Q)
        num, den = c // g, Q // g
        for q in continued_fraction_denominators(num, den):
            if q < 2:
                continue
            for m in [1, 2, 3, 4, 8, 16]:
                if q * m == r:
                    return True, r, (c, q, m)
    return False, r, None


def find_a(N, Q=256):
    """
    Return list of all room-compatible (a, r, c, q, m) tuples for factoring N.
    Returns empty list if none found.
    """
    result = []
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
        if math.gcd(a, N) != 1:
            continue
        ok, r, detail = room_compatible(N, a, Q)
        if ok:
            c, q, m = detail
            result.append((a, r, c, q, m))
    return result


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    N = int(sys.argv[1])
    banks = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 1
    show_all = "--all" in sys.argv

    D = 256
    Q = banks * D

    if N < 4:
        print(f"N={N} is too small to factor meaningfully.")
        sys.exit(1)

    # Quick GCD check for trivial factors
    for a in [2, 3, 5, 7, 11]:
        g = math.gcd(a, N)
        if 1 < g < N:
            print(f"N={N} shares factor {g} — trivially factored: {g} × {N//g}")
            sys.exit(0)

    if show_all:
        pairs = find_a(N, Q)
        print(f"All room-compatible (N={N}, a) pairs with Q={Q} (banks={banks}):")
        print(f"{'a':>4s} {'r':>6s} {'peak':>8s} {'q':>6s} {'×m':>4s} {'=r':>6s}")
        print("─" * 44)
        if not pairs:
            print("  (none — try more banks)")
        for a, r, c, q, m in pairs:
            print(f"{a:4d} {r:6d} {c:8d} {q:6d} {m:4d} {r:6d}")
    else:
        result = find_a(N, Q)
        if not result:
            # Try with more banks as fallback
            for extra_banks in [2, 4, 8, 16, 32, 64, 128]:
                Q2 = extra_banks * D
                result = find_a(N, Q2)
                if result:
                    a, r, c, q, m = result[0]
                    print(f"N={N}: a={a}  r={r}  (needs banks={extra_banks}, Q={Q2})")
                    print(f"  peak at bin {c} → {c}/{Q2} cf-denom q={q} ×{m}={r}")
                    sys.exit(0)
            print(f"N={N}: no compatible a found (try more banks or different Q)")
            sys.exit(1)
        else:
            a, r, c, q, m = result[0]
            g = math.gcd(c, Q)
            print(f"N={N}  a={a}  r={r}  Q={Q}  banks={banks}")
            print(f"  QFT peak: bin {c} ({c}/{Q} = {c//g}/{Q//g})")
            print(f"  Continued fraction: q={q} ×{m} = {q*m} (=r)")


if __name__ == "__main__":
    main()
