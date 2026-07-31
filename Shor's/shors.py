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
    # Match PROB output: |k>=prob or |k⟩=prob
    ps = {}
    for m in re.finditer(r'\|(\d+)[⟩>]=(\d+\.\d+)', output):
        ps[int(m.group(1))] = float(m.group(2))
    return ps

def griffiths_niu(N, a, D, freq, rate, gain, trials, use_sdr):
    M = 2*N
    print(f"[gn] N={N} a={a} D={D} {'room' if use_sdr else 'soft'}")
    
    sqrt_n = int(math.isqrt(N))
    remainder = N - sqrt_n * (N // sqrt_n)
    max_k = 2 * N.bit_length() + 6
    
    for trial in range(trials):
        phase = 0.0
        handoff_x = sqrt_n
        prev_phi = -1
        stable_count = 0
        
        # Per-iteration subprocess (clean WF each time), early exit on geodesic/phase
        for k in range(max_k):
            c = mod_pow(a, 1 << k, N)
            if c == 0 or gcd(c, N) != 1:
                phase = phase / 2.0; continue
            angle = (2.0 * math.pi * c / N) % (2.0 * math.pi)
            script = [
                "RESET", "SET 2 0.70710678", "SET 3 0.70710678",
                f"Z {angle:.15f} 3", "HCTRL", "PROB", "QUIT"
            ]
            out = run_engine(script, D, freq, rate, gain)
            if not out: print(f"  k={k} fail"); break
            ps = parse_probs(out)
            p0 = ps.get(2, 0); p1 = ps.get(3, 0)
            tot = p0 + p1
            if tot < 1e-12: p0 = p1 = 0.5; tot = 1.0
            p0 /= tot; p1 /= tot
            
            phi_k = 2.0 * math.acos(max(-1.0, min(1.0, math.sqrt(p0))))
            phase = (phase + phi_k) / 2.0
            
            # Geodesic
            z_val = int(phi_k * N / (2.0 * math.pi))
            resonance = (c * handoff_x - z_val * remainder) % N
            if resonance == 0: resonance = c % N
            g = gcd(resonance, N)
            if 1 < g < N and N % g == 0:
                print(f"  GEODESIC k={k}: c={c} φ={phi_k/math.pi:.4f}π → gcd={g}")
                print(f"[gn] ★ {N} = {min(g,N//g)} × {max(g,N//g)}")
                return (min(g,N//g), max(g,N//g))
            handoff_x = (handoff_x * c + z_val) % N
            
            # Early exit on phase convergence
            if prev_phi > 0 and abs(phi_k - prev_phi) < 0.001:
                stable_count += 1
            else:
                stable_count = 0
            prev_phi = phi_k
            
            if k < 12 or k % 8 == 0:
                print(f"  k={k:3d} c={str(c)[:20]}... p0={p0:.4f} p1={p1:.4f} "
                      f"φ_k={phi_k/math.pi:.4f}π phase={phase/math.pi:.6f}π")
            
            if k >= N.bit_length() and phase > 1e-8:
                # Extract period via continued fractions on accumulated phase
                frac = (phase / (2.0 * math.pi)) % 1.0
                if frac > 1e-10:
                    p0, q0, p1, q1 = 0, 1, 1, 0
                    x = frac
                    for _ in range(60):
                        if abs(x) < 1e-15: break
                        ai = int(math.floor(x + 1e-12))
                        if ai < 0: break
                        p2, q2 = ai * p1 + p0, ai * q1 + q0
                        if q2 > N or q2 < 0: break
                        if q2 > 1 and mod_pow(a, q2, N) == 1:
                            r = q2
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
                        p0, q0, p1, q1 = p1, q1, p2, q2
                        x = 1.0 / (x - ai)
                print(f"  k={k:3d} c={str(c)[:20]}... p0={p0:.4f} p1={p1:.4f} "
                      f"φ_k={phi_k/math.pi:.4f}π phase={phase/math.pi:.6f}π")
            
            # Phase convergence → try period
            if stable_count >= 10:
                print(f"  Phase stable at {phase/math.pi:.6f}π after {k} iterations, early exit")
                break
            
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
