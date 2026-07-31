#!/usr/bin/env python3
"""shors.py — Griffiths-Niu via CR_BIN phase rotations + SDR EXP."""

import subprocess, sys, math, random, re, os, tempfile, argparse
from math import gcd

ENGINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sdr_ether")

def mod_pow(b,e,m):
    r=1;b%=m
    while e:
        if e&1: r=(r*b)%m
        b=(b*b)%m; e>>=1
    return r

def pick_coprime(N,avoid=None):
    avoid=avoid or set()
    for _ in range(100):
        a=random.randint(2,N-2)
        if a not in avoid and gcd(a,N)==1: return a
    return None

def run_engine(script, D, freq, rate, gain):
    s='\n'.join(script)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.qvm', delete=False) as f:
        f.write(s); sp=f.name
    try:
        r=subprocess.run([ENGINE,str(D),str(freq),str(rate),str(gain),'--vm',sp],
                        capture_output=True,text=True,timeout=300,env={**os.environ})
        return r.stdout+r.stderr
    except: return None
    finally: os.unlink(sp)

def parse_probs(output):
    ps = {}
    for m in re.finditer(r'\|(\d+)[⟩>]=(\d+\.\d+)', output):
        ps[int(m.group(1))] = float(m.group(2))
    return ps
def extract_period_int(P, K, N, a):
    """Extract period r from exact rational P/2^K via integer continued fractions."""
    if K < 1 or P == 0:
        return None
    # Integer Euclidean algorithm on (num, den) — no float conversions
    num, den = P, 1 << K
    n0, d0, n1, d1 = 0, 1, 1, 0
    while den > 0:
        ai = num // den
        if ai < 0: break
        n2, d2 = ai * n1 + n0, ai * d1 + d0
        if d2 > N or d2 < 0: break
        if d2 > 1 and mod_pow(a, d2, N) == 1:
            return d2
        n0, d0, n1, d1 = n1, d1, n2, d2
        num, den = den, num - ai * den
        if den == 0: break
    return None

def griffiths_niu(N, a, D, freq, rate, gain, trials, use_sdr):
    M = 2*N
    print(f"[gn] N={N} a={a} D={D} {'room' if use_sdr else 'soft'}")
    
    sqrt_n = int(math.isqrt(N))
    remainder = N - sqrt_n * (N // sqrt_n)
    max_k = 2 * N.bit_length() + 6
    
    for trial in range(trials):
        phase = 0.0
        # Integer-phase accumulator: exact binary fraction
        P = 0  # numerator, denominator = 2^k
        K = 0  # number of bits accumulated
        handoff_x = sqrt_n
        
        # Per-iteration subprocess (clean WF each time), early exit on geodesic/phase
        # Per-iteration: feedback depends on prior measurement outcomes
        for k in range(max_k):
            c = mod_pow(a, 1 << k, N)
            if c == 0 or gcd(c, N) != 1:
                P <<= 1; K += 1; phase = phase / 2.0; continue
            angle = (2.0 * math.pi * c / N) % (2.0 * math.pi)
            script = ["RESET", "SET 2 0.70710678", "SET 3 0.70710678"]
            if K > 0:
                fb = (math.pi * P) / (1 << K)
                script.append(f"Z {-fb/2:.15f} 2")
                script.append(f"Z {+fb/2:.15f} 3")
            script.extend([f"Z {angle:.15f} 3", "HCTRL", "PROB", "QUIT"])
            out = run_engine(script, D, freq, rate, gain)
            if not out: print(f"  k={k} fail"); break
            ps = parse_probs(out)
            p0 = ps.get(2, 0); p1 = ps.get(3, 0)
            tot = p0 + p1
            if tot < 1e-12: p0 = p1 = 0.5; tot = 1.0
            p0 /= tot; p1 /= tot
            
            phi_k = 2.0 * math.acos(max(-1.0, min(1.0, math.sqrt(p0))))
            phase = (phase + phi_k) / 2.0
            
            # Integer-phase: store bit with full precision (binary threshold)
            bit = 1 if phi_k > math.pi/2 else 0
            P = (bit << K) | P
            K += 1
            
            # Geodesic — pad z_val ±1 in SDR mode for phase jitter
            z_val = int(phi_k * N / (2.0 * math.pi))
            for dz in (range(-1, 2) if use_sdr else [0]):
                zv = z_val + dz
                resonance = (c * handoff_x - zv * remainder) % N
                if resonance == 0: resonance = c % N
                g = gcd(resonance, N)
                if 1 < g < N and N % g == 0:
                    print(f"  GEODESIC k={k}: c={c} φ={phi_k/math.pi:.4f}π → gcd={g}")
                    print(f"[gn] ★ {N} = {min(g,N//g)} × {max(g,N//g)}")
                    return (min(g,N//g), max(g,N//g))
            handoff_x = (handoff_x * c + z_val) % N
            
            if k < 12 or k % 8 == 0:
                print(f"  k={k:3d} c={str(c)[:20]}... p0={p0:.4f} p1={p1:.4f} "
                      f"φ_k={phi_k/math.pi:.4f}π phase={phase/math.pi:.6f}π")
            
            # Period extraction via continued fractions on integer phase
            if k >= N.bit_length() and K > 0:
                r = extract_period_int(P, K, N, a)
                if r:
                    print(f"  Converged: r={r} (k={k})")
                    if r % 2 == 0:
                        h = mod_pow(a, r // 2, N)
                        if h != N - 1:
                            f1, f2 = gcd(h+1, N), gcd(h-1, N)
                            if 1 < f1 < N and 1 < f2 < N and f1*f2 == N:
                                print(f"[gn] ★ {N} = {min(f1,f2)} × {max(f1,f2)} (phase r={r})")
                                return (min(f1,f2), max(f1,f2))
                    print(f"  r={r} no factor")
                    return None
        # Post-loop
        if K > 0:
            r = extract_period_int(P, K, N, a)
            if r:
                print(f"  Post-loop r={r}")
                if r % 2 == 0:
                    h = mod_pow(a, r // 2, N)
                    if h != N - 1:
                        f1, f2 = gcd(h+1, N), gcd(h-1, N)
                        if 1 < f1 < N and 1 < f2 < N and f1*f2 == N:
                            print(f"[gn] ★ {N} = {min(f1,f2)} × {max(f1,f2)} (post-loop)")
                            return (min(f1,f2), max(f1,f2))
                print(f"  r={r} no factor")
                return None
        print(f"  Trial complete: phase={phase/math.pi:.6f}π (no period found)")
    return None

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('N', type=int, nargs='?', default=0)
    p.add_argument('a', type=int, nargs='?', default=None)
    p.add_argument('--D', type=int, default=32768)
    p.add_argument('--freq', type=float, default=150e6)
    p.add_argument('--rate', type=float, default=2048000)
    p.add_argument('--gain', type=int, default=496)
    p.add_argument('--trials', type=int, default=3)
    p.add_argument('--sdr', action='store_true', default=True)
    p.add_argument('--soft', action='store_true')
    args = p.parse_args()
    if args.N < 2: sys.exit("N required")
    use_sdr = args.sdr and not args.soft
    
    if not use_sdr:
        subprocess.run(['sudo','rmmod','rtl2832_sdr'], capture_output=True, timeout=5)
    try:
        for t in range(args.trials):
            if args.a: av, args.a = args.a, None
            else: av = pick_coprime(args.N)
            if not av: break
            if gcd(av, args.N) != 1: print(f"gcd={gcd(av,args.N)}"); break
            print(f"[gn] T{t+1}: a={av}")
            r = griffiths_niu(args.N, av, args.D, args.freq, args.rate, args.gain, args.trials, use_sdr)
            if r: break
    finally:
        if not use_sdr:
            subprocess.run(['sudo','modprobe','rtl2832_sdr'], capture_output=True, timeout=5)
