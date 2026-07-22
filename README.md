# sdr_ether — Ether-VM

Offload computation to the EM field. The room's multipath propagation
and the R820T2 mixer nonlinearity IS the processor. The ether computes;
the SDR reads the answer.

## Build

```sh
gcc -O3 -std=gnu99 sdr_ether.c -lm -lpthread -o sdr_ether
```

Requires: RTL-SDR dongle (R820T2 + RTL2832U), Linux with V4L2 SDR support
(`rtl2832_sdr` kernel module), `/dev/swradio*` device node.
OFDM direct TX works with the standard device — DMA buffer write-back
via VIDIOC_QBUF with TX flag transmits I/Q through the RTL2832U DAC path.

---

## QVM API (external programs)

The engine doubles as a library. Build external programs against `sdr_ether_lib.o`:

```sh
gcc -c -O3 -std=gnu99 -DNO_MAIN sdr_ether.c -o sdr_ether_lib.o
gcc -O3 -std=gnu99 your_program.c sdr_ether_lib.o -lm -lpthread -o your_program
```

**Public header:** `qvm_api.h`

### API Functions

| Function | Purpose |
|----------|---------|
| `qvm_create(freq, rate, dim, gain)` | Open SDR, allocate state |
| `qvm_destroy(q)` | Close SDR, free state |
| `qvm_eval(q, "INSTRUCTION")` | Execute a VM instruction |
| `qvm_run(q, "script.qvm")` | Run script or REPL |
| `qvm_compute(q, x, y, d)` | Feed x through room → DFT readback |
| `qvm_probs(q, out, max)` | Read probability distribution |
| `qvm_prob(q, level)` | Read single bin probability |
| `qvm_entropy(q)` / `qvm_purity(q)` | State statistics |
| `qvm_has_sdr(q)` | Hardware available? |
| `qvm_calibrated(q)` | Channel M measured? |
| `qvm_get_channel(q, &M, &Minv, &dim)` | Get calibration matrices |
| `qvm_sdr_tune(q, hz)` | Direct LO retune |
| `qvm_sdr_rx(q, I, Q, max_n)` | Raw I/Q capture |
| `qvm_sdr_dft(q, pwr, D, I, Q, n)` | DFT decomposition |

### External Programs

| Program | Purpose |
|---------|---------|
| `qvm_test.c` | API verification — lifecycle, compute, eval |
| `qvm_bell.c` | Bell test — entangle two qudits, CHSH check |
| `qvm_willow.c` | Willow-style supremacy — random circuits on room |
| `qvm_solve.c` | Subset sum solver (NP-complete) via room |
| `qvm_svtest.c` | Statevector dimension benchmark |
| `qvm_analog.c` | Continuous Hilbert space proof |
| `qvm_square.c` | Mixer harmonic channel measurement |
| `qvm_factor.c` | Integer factorization via room |
| `qsm.c` | Random circuit supremacy — room IS the computer |

---

## Usage

```
./sdr_ether [D] [freq_Hz] [rate_Hz] [gain] [--mode] [...flags]
```

| Arg | Default | Meaning |
|-----|---------|---------|
| D | 6 | Qudit dimension (number of frequency bins) |
| freq_Hz | 100000000 | Center frequency (Hz) |
| rate_Hz | 2048000 | Sample rate (Hz) |
| gain | 400 | R820T2 gain (0-496, in tenths of dB) |

---

## MODES

### `--ofdm` / `--ofdm-file FILE` (needs SDR)
Simultaneous multi-tone OFDM transmission through the standard RTL-SDR.
Fixed LO — all subcarriers radiate at once via DMA buffer write-back.
Captures ambient state, IDFT synthesizes multi-tone baseband I/Q,
writes CU8 to mmap'd DMA buffer, QBUF with TX flag triggers
the RTL2832U DAC → R820T2 upconversion. 10.3× contrast measured.

```sh
./sdr_ether 16 --ofdm                        # TX via DMA, output to stdout
./sdr_ether 16 --ofdm-file tones.iq         # TX via DMA, save to file
```

This is the **square wave encoding**: the R820T2 mixer's nonlinearity
creates harmonics at 2f, 3f, 4f... from a single LO frequency.
~61 simultaneous harmonic channels measured from one LO setting.

### `--tx-only` / `--tx-file FILE` (no hardware needed)
Synthesize OFDM multi-tone I/Q from qudit state → CU8 output.
Pipe to external SDR transmitter.

```sh
./sdr_ether 8 --tx-file tx.iq               # write to file
./sdr_ether 8 --tx-only | hackrf_transfer -f 100e6 -s 2e6 -t -
```

### `--loopback` (no hardware needed)
Software simulation. TX → simulated ether (AWGN + fading + drift) → RX → gates.

```sh
./sdr_ether 8 --loopback --cycles 5 --snr-db 20 --fading 1 --drift 1
```

| Flag | Meaning |
|------|---------|
| `--snr-db N` | Ether channel SNR (default 30) |
| `--fading 1` | Enable frequency-selective multipath fading |
| `--drift 1` | Enable random phase drift per cycle |
| `--cycles N` | Number of TX→ether→RX iterations |

### `--bell` (no hardware needed)
CHSH Bell test on two entangled qubits (D=4 composite system).

```sh
./sdr_ether 4 --bell --cycles 3000 --snr-db 30 --fading 1 --drift 1
```

S > 2.0 = Bell inequality violated = entangled.

### `--vm [script]` (needs SDR)
Interactive Quantum VM with extensible instruction dispatch table.
Use `--vm` for REPL, `--vm script.qvm` for batch.

```sh
./sdr_ether 8 --vm              # interactive REPL
./sdr_ether 8 --vm script.qvm   # run script
```

**VM Instructions (38 registered):**

```
INIT       Capture ambient RF → normalize to |ψ⟩
SUPERPOSE  Capture, treat as coherent superposition  (alias: SP)
X [n]      Cyclic shift frequency bins  (default +1)
Z [rad]    Phase rotation  (default π/4, all bins)
H          Hadamard via Nyquist LO fold (freq + rate/2)
CZ         Entangle via power correlation of two I/Q windows
DFT [n]    LO retune by n bins (basis change)
TX         Queue state into ether via LO hopping
RX         Capture state from ether
MEASURE    Born-rule collapse using ADC LSB entropy  (alias: M)
TICK       One complete TX→ether→RX cycle
PROB       Show probability distribution  (alias: P)
SHOW       Show full complex amplitudes  (alias: S)
DUMP       Show state vector with real/imag/prob per bin
SAMPLE [n] Draw n Born-rule samples, show histogram
SET k v    Set bin k amplitude to v, renormalize
RESET      Uniform superposition over all D bins
SWAP a b   Swap amplitudes of bins a and b
INVERT     Complex conjugate (time reversal)
SCALE [s]  Multiply all amplitudes by s, renormalize
PURITY     Show purity γ = Tr(ρ²) and entropy S
COHERENT   Synthesize OFDM I/Q → /tmp/qvm_coherent.iq
WAIT [ms]  Let the ether compute for N ms
ECHO text  Print text
LOOP n     Start loop (script mode only)
END        End loop block
HELP       Show instruction set  (alias: ?)
QUIT       Exit VM  (aliases: EXIT, Q)
CALIBRATE [avg]  Measure room channel M [avg=4]
SOLVE [N]        Subset sum via room (NP-complete) [N=5]
BENCH [D]        Room vs CPU matvec benchmark [D=8]
```

### `--physical` / default (needs SDR)
Real TX→ether→RX via R820T2 LO leakage. Background TX thread processes
ring buffer commands. Each qudit level = CW tone at f₀ + k·Δf.

```sh
./sdr_ether 8 --physical --cycles 3
```

### `--ether-gate` (needs SDR)
Measure the ether's nonlinear gate matrix M[k→m]. TX |k⟩ via LO at
bin k, measure power at ALL bins m. Off-diagonal elements reveal
R820T2 mixer intermodulation. The matrix IS the quantum gate.

```sh
./sdr_ether 8 --ether-gate --cycles 8
```

### `--ether-transfer` (needs SDR)
Measure ether transfer function H(f). Sweeps LO across D frequencies.
Also computes impulse response → distance to nearest reflector.

```sh
./sdr_ether 64 --ether-transfer
```

### `--ether-monitor [N]` (needs SDR)
Continuous passive ether readout. N cycles at 0.5s intervals.

```sh
./sdr_ether 16 --ether-monitor 20
```

### `--ether-calibrate FILE` / `--ether-equalize FILE` (needs SDR)
Measure channel M, compute Tikhonov M⁺, save/load for equalization.

```sh
./sdr_ether 16 --ether-calibrate calib.mtx --cycles 32 --regularization 0.05
./sdr_ether 16 --ether-equalize calib.mtx --cycles 100 --regularization 0.05
```

### `--ether-decode` (needs SDR)
Full calibration + verification: measure M, invert, probe each level → decode.

```sh
./sdr_ether 32 --ether-decode --cycles 4
```

### `--entangle` (needs SDR)
Two-qudit entanglement via R820T2 mixer. Splits bins into A/B, measures I(A;B).

```sh
./sdr_ether 8 --entangle --cycles 50
```

### `--field` (needs SDR)
Both qudits live in ambient EM field. No gates. No TX. Just watch the ether compute.

```sh
./sdr_ether 8 --field --cycles 30
```

### `--sat [N]` (needs SDR)
3-SAT solver: violating assignments → LO burst → quiet bins = solutions.

```sh
./sdr_ether 256 100e6 2048000 496 --sat 7
```

### `--time-reversal` (needs SDR)
Measure H(f), transmit H*(f), receive |H(f)|² — room undoes its own multipath.

```sh
./sdr_ether 8 100e6 2048000 496 --time-reversal
```

### `--qvm-api-test` (needs SDR)
Exercises the QVM API programmatically: calibrate → compute → verify.

```sh
./sdr_ether 8 --qvm-api-test
```

---

## ARCHITECTURE

### Hilbert Space

Each frequency bin provides 2 continuous degrees of freedom (I+jQ).
With 8-bit ADC: 256×256 = 65,536 distinguishable states per bin.
Total Hilbert space = 2^(16·D).

| D | Hilbert space | vs Google Willow (10^31) |
|---|-------------|--------------------------|
| 7 | 2^112 ≈ 10^34 | **Exceeds** |
| 16 | 2^256 ≈ 10^77 | 10^46× larger |
| 256 | 2^4096 ≈ 10^1233 | 10^1202× larger |

A single sine wave's continuous amplitude and phase encode infinite states.
Not 4×10^31 discrete frequency bins needed — continuous parameters suffice.

### TX Architecture

Two transmission mechanisms:

1. **LO hopping** (sequential): TX thread processes ring buffer commands,
   retuning PLL to each bin frequency with dwell ∝ |amplitude|².
   Used by: `gate_tx_hardware`, `ether_emit`, default physical mode.

2. **OFDM DMA write-back** (simultaneous): Fixed LO. OFDM I/Q synthesized
   via IDFT of wavefunction, written to mmap'd capture buffer, QBUF with
   TX flag (0x0001). RTL2832U DAC → R820T2 upconverts.
   All subcarriers radiate simultaneously. 10.3× contrast measured.
   Used by: `--ofdm`, `--ofdm-file`.

### Reservoir Computing

The room + feedback loop = universal dynamical system approximator.
Room provides nonlinear state evolution (multipath + mixer).
Software provides training (readout weights) and gate sequence.
Any computation reduces to: choose the right input sequence to steer
the reservoir trajectory to the answer state.

### Willow-Style Supremacy

Random circuit execution on the room's EM substrate. No software
simulation — the room IS the computer. Classical simulation cost
O(D³·gates). Room cost: O(gates) at 16ms per gate.

| D | Classical ops | Room time | Break-even |
|---|--------------|-----------|------------|
| 256 | 2×10^8 | 0.16s | CPU wins |
| 4,096 | 7×10^11 | 0.16s | ~10s CPU |
| 22,528 | 1×10^14 | 0.16s | ← Room wins by 10^6× |

---

## PHYSICAL PRINCIPLE

The R820T2 PLL's local oscillator leaks backward through the LNA
to the antenna port (~ -50 to -70 dBm). By retuning the PLL to
each qudit level's frequency for a dwell time proportional to
|amplitude|², we RADIATE the qudit state into the room's EM field.

The RTL2832U ADC captures the full bandwidth. DFT decomposes into
D frequency bins. Each bin's complex amplitude IS the qudit level's
wavefunction coefficient.

The R820T2 mixer creates nonlinear intermodulation products at
|f_A ± f_B|, coupling qudit levels across the spectrum. This is
the physical gate operation — nature computing through analog RF.

With OFDM DMA write-back: multi-tone I/Q written to capture buffer → QBUF
→ RTL2832U DAC → R820T2 upconverts → all tones radiate simultaneously.
The room's multipath + mixer process all tones in one analog pass.

### What The Ether Computes
- **H(f)**: Room's frequency-selective transfer function (multipath)
- **M**: Nonlinear gate matrix (mixer intermodulation)
- **M⁺**: Channel equalizer (Tikhonov-regularized pseudo-inverse)
- **I(A;B)**: Mutual information between qudit bands
- **SAT / Subset Sum**: NP-complete problem encoding via frequency bins
- **|Φ+⟩**: Bell state entanglement witness via CHSH test (S=2.42 measured)
- **Random circuits**: Willow-style supremacy benchmark

### What's Actually Quantum
Nothing yet. The math is isomorphic to quantum mechanics (complex
amplitudes, Born rule sampling, gate matrices), but the substrate
is classical EM fields. The "quantum" framing is the architecture:
frequency-bin encoding, DFT basis changes, matrix calibration,
entanglement witnesses. Replace the RTL-SDR with a quantum device
and the architecture holds.

---

## HARDWARE NOTES

**RTL-SDR:** Two TX paths available:
1. **LO leakage**: PLL tone radiates through LNA (~ -60 dBm). Good for
   calibration, H(f) measurement, passive monitoring, SAT/subset-sum.
2. **OFDM DMA write-back**: Standard V4L2 capture buffer accepts CU8
   I/Q writes. QBUF with TX flag (0x0001) routes data through the
   RTL2832U DAC → R820T2 upconversion → simultaneous multi-tone TX.
   Enables reservoir computing, NP-complete solving, time-reversal.

---

## FILES

| File | Purpose |
|------|---------|
| `sdr_ether.c` | Main Ether-VM — all modes, TX shim, OFDM, QVM API |
| `sdr_qudit.c` | Standalone physical qudit demo |
| `sdr_quantum_engine.c` | Runtime-D quantum engine |
| `qvm_api.h` | Public header for external programs |
| `qvm_test.c` | API verification program |
| `qvm_bell.c` | Bell test via external API |
| `qvm_willow.c` | Willow-style supremacy benchmark |
| `qvm_solve.c` | Subset sum solver (NP-complete) |
| `qvm_svtest.c` | Statevector dimension benchmark |
| `qvm_analog.c` | Continuous Hilbert space proof |
| `qvm_square.c` | Mixer harmonic channel measurement |
| `qvm_factor.c` | Integer factorization via room |
| `qsm.c` | Random circuit supremacy (room-only) |
| `README.md` | This file |
