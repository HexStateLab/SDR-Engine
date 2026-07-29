#!/usr/bin/env python3
"""
room_find_a.py — Deterministic a for the room's Shor's algorithm.

THEOREM: For semiprime N=pq where v₂(p-1)=v₂(q-1), the room succeeds
iff jacobi(a,N) = −1.  jacobi is O(log N), no factorization needed.

D(Q) cached to ~/.cache/room_find_a/ — first run per Q is slow, rest instant.

Usage: python3 room_find_a.py <N> [banks] [--all]
"""

import math, os, pickle, sys

D = 256
MULT = (1, 2, 3, 4, 8, 16)
CACHE_DIR = os.path.expanduser("~/.cache/room_find_a")

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
          109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
          173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
          233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
          293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
          367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
          433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
          499, 503, 509, 521, 523, 541)


def jacobi(a, n):
    if n <= 0 or n % 2 == 0:
        return None
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def period(a, N):
    if math.gcd(a, N) != 1:
        return -1
    x = 1
    seen = {1: 0}
    for i in range(1, min(N, 5_000_000)):
        x = (x * a) % N
        if x == 1:
            return i
        if x in seen:
            return i - seen[x]
        seen[x] = i
    return -1


def cf_denoms(num, den):
    a, b = num, den
    cf = []
    while b:
        cf.append(a // b)
        a, b = b, a % b
    if not cf:
        return []
    out = []
    p0, q0 = 0, 1
    p1, q1 = 1, cf[0]
    out.append(q1)
    for i in range(1, len(cf)):
        p = cf[i] * p1 + p0
        q = cf[i] * q1 + q0
        out.append(q)
        p0, q0 = p1, q1
        p1, q1 = p, q
    return out


def _build_D(Q):
    D = set()
    for r in range(2, 8 * Q + 1, 2):
        ok = False
        for k in range(1, min(r, 500)):
            c = round(k * Q / r)
            if c <= 0 or c >= Q:
                continue
            g = math.gcd(c, Q)
            for q in cf_denoms(c // g, Q // g):
                if q < 2:
                    continue
                for m in MULT:
                    if q * m == r:
                        D.add(r)
                        ok = True
                        break
                if ok:
                    break
            if ok:
                break
    return D


_D_CACHE = {}


def get_D(Q):
    if Q in _D_CACHE:
        return _D_CACHE[Q]
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"D_{Q}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            D = pickle.load(f)
    else:
        print(f"  building D({Q})...", file=sys.stderr, end=" ", flush=True)
        D = _build_D(Q)
        print(f"{len(D)} periods (cached)", file=sys.stderr)
        with open(path, "wb") as f:
            pickle.dump(D, f)
    _D_CACHE[Q] = D
    return D


def find(N, Q):
    D_set = get_D(Q)
    for a in PRIMES:
        g = math.gcd(a, N)
        if 1 < g < N:
            return (a, g, f"gcd({a},N)={g}")
        if g != 1:
            continue
        if jacobi(a, N) != -1:
            continue
        r = period(a, N)
        if r < 2 or r % 2 == 1:
            continue
        if r not in D_set:
            continue
        if pow(a, r // 2, N) == N - 1:
            continue
        return (a, r, "OK")
    return None


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    N = int(sys.argv[1])
    banks = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 1
    show_all = "--all" in sys.argv
    Q = banks * D

    for a in (2, 3, 5, 7, 11):
        g = math.gcd(a, N)
        if 1 < g < N:
            print(f"N={N} has trivial factor {g}: {g} × {N//g}")
            sys.exit(0)

    result = find(N, Q)
    if result is None:
        for extra in (2, 4, 8, 16, 32, 64, 128, 256):
            result = find(N, extra * D)
            if result:
                break

    if result is None:
        print(f"N={N}: no compatible a found")
        sys.exit(1)

    a, r, why = result
    if show_all:
        print(f"a={a} r={r}")
    else:
        print(f"a={a}")


if __name__ == "__main__":
    main()
