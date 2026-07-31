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
    m=re.search(r'\[PROB\](.*?)S=',output)
    if not m: return {}
    ps={}
    for p in m.group(1).split():
        kv=re.match(r'\|(\d+)[⟩>]=(\d+\.\d+)',p)
        if kv: ps[int(kv.group(1))]=float(kv.group(2))
    return ps

def griffiths_niu(N, a, D, freq, rate, gain, trials, use_sdr):
    """Griffiths-Niu via CR_BIN phase rotations on Fourier bins."""
    
    M = 2*N  # unused now — CPHASE/HCTRL/JMEAS scale to any N
    print(f"[gn] N={N} a={a} D={D} {'room' if use_sdr else 'soft'}")
    
    for trial in range(trials):
        phase = 0.0
        # Geodesic accumulator state (from Squarer.py)
        sqrt_n = int(math.isqrt(N))
        remainder = N - sqrt_n * (N // sqrt_n)
        handoff_x = sqrt_n  # initial coordinate shadow
        
        for k in range(2 * N.bit_length() + 6):
            c = mod_pow(a, 1 << k, N)
            if c == 0 or gcd(c, N) != 1:
                phase = phase / 2.0; continue
            
            # Circuit: QFT → CR_BIN phase rotations → IQFT → H on ctrl → measure
            # Circuit: CPHASE → HCTRL → JMEAS (no QFT needed!)
            # CPHASE applies e^{i·2π·c·x/N} to |1⟩|x⟩ components
            # This is the eigenvalue directly — no Fourier basis needed
            script = [
                "RESET",
                "SET 2 0.70710678",  # |0⟩|1⟩  (joint-state bins 2,3)
                "SET 3 0.70710678",  # |1⟩|1⟩
                f"CPHASE {c} {N}",   # controlled phase = eigenvalue
                f"HCTRL {N}",        # H on control → interference
                f"JMEAS {N}",        # measure
                "QUIT"
            ]
            
            out = run_engine(script, D, freq, rate, gain)
            if not out: print(f"  k={k} fail"); break
            
            m = re.search(r'JMEAS.*P\(.*?\)=(\d+\.\d+).*?P\(.*?\)=(\d+\.\d+)', out)
            if not m: print(f"  k={k} no JMEAS"); continue
            p0, p1 = float(m.group(1)), float(m.group(2))
            
            phi_k = 2.0 * math.acos(max(-1.0, min(1.0, math.sqrt(p0))))
            phase = (phase + phi_k) / 2.0
            
            # Geodesic resonance check — direct GCD extraction
            z_val = int(phi_k * N / (2.0 * math.pi))  # phase → z-coordinate
            resonance = abs(c * handoff_x - z_val * remainder)
            g = gcd(resonance, N)
            if 1 < g < N and N % g == 0:
                f1, f2 = min(g, N//g), max(g, N//g)
                print(f"  GEODESIC: c={c} φ={phi_k/math.pi:.4f}π → gcd={g}")
                print(f"[gn] ★ {N} = {f1} × {f2}")
                return (f1, f2)
            # Accumulate coordinate shadow
            handoff_x = (handoff_x * c + z_val) % N
            
            if k < 12 or k % 8 == 0:
                print(f"  k={k:3d} c={c:5d} p0={p0:.4f} p1={p1:.4f} "
                      f"φ_k={phi_k/math.pi:.4f}π phase={phase/math.pi:.6f}π")
            
            if phase > 1e-8:
                r = round(1.0 / (phase / (2.0 * math.pi)))
                if r >= 2 and r <= N and mod_pow(a, r, N) == 1:
                    print(f"  Found: r={r}")
                    if r % 2 == 0:
                        h = mod_pow(a, r // 2, N)
                        if h != N - 1:
                            f1, f2 = gcd(h + 1, N), gcd(h - 1, N)
                            if 1 < f1 < N and 1 < f2 < N and f1 * f2 == N:
                                print(f"[gn] ★ {N} = {min(f1,f2)} × {max(f1,f2)}")
                                return (min(f1, f2), max(f1, f2))
                    return None
        
        # End of k-loop: try accumulated phase
        if phase > 1e-8:
            r = round(1.0 / (phase / (2.0 * math.pi)))
            if r >= 2 and r <= N and mod_pow(a, r, N) == 1:
                if r % 2 == 0:
                    h = mod_pow(a, r // 2, N)
                    if h != N - 1:
                        f1, f2 = gcd(h + 1, N), gcd(h - 1, N)
                        if 1 < f1 < N and 1 < f2 < N and f1 * f2 == N:
                            print(f"[gn] ★ {N} = {min(f1,f2)} × {max(f1,f2)} (final phase)")
                            return (min(f1, f2), max(f1, f2))
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
