# ether_boost — R820T2 TX Register Boost for Non-TX-Native SDRs

## Why This Exists

Standard RTL-SDR dongles with the R820T2 tuner are **receive-only** by design.
The signal path flows LNA → mixer → IF, and all internal registers are configured
for low-noise reception with minimal LO leakage.

The Ether-VM (`sdr_ether`) transmits OFDM subcarriers via **DMA buffer write-back**:
CU8 I/Q samples are written to mmap'd V4L2 capture buffers and submitted with
QBUF + TX flag (0x0001). The RTL2832U routes this data through its internal DAC,
mixing the baseband I/Q up to the R820T2's LO frequency.

For this reverse signal path to work at usable power, the R820T2 must be
reconfigured to: bypass the LNA (reverse gate), maximize mixer bias current
for bidirectional conduction, and drive the VCO at maximum amplitude so the
LO leaks backward through the LNA to the antenna port.

**Without this boost, TX power is ~-70 dBm — too weak for the room's multipath.**
**With the boost, TX power increases ~15-20 dB, enough for measurable intermodulation.**

## Files

| File | Purpose |
|------|---------|
| `ether_boost.c` | User-space R820T2 register poke via USB control transfer |
| `boost.sh` | Auto-detect device, set permissions, build & run ether_boost |

## Usage

### Option 1: ether_boost (recommended, no kernel build)

```sh
# One-time: find and unlock the USB device
lsusb -d 0bda:2838
sudo chmod 666 /dev/bus/usb/001/016   # adjust bus/device numbers

# Build and run
gcc -O2 ether_boost.c -o ether_boost
./ether_boost
```

### Option 2: boost.sh (auto-detect)

```sh
./boost.sh
```

Finds the RTL2838 automatically, sets USB permissions, builds and runs ether_boost.

## What It Does

Pokes 5 R820T2 registers via I2C (address 0x34) over USB vendor control transfer:

| Register | Value | Effect |
|----------|-------|--------|
| 0x05 | 0x00 | LNA gain bypass (disables forward gain, enables reverse path) |
| 0x06 | 0x01 | LNA bypass enable |
| 0x07 | 0x1F | Maximum mixer bias current (bidirectional conduction) |
| 0x0A | 0x60 | Maximum VCO amplitude |
| 0x10 | 0x0F | Maximum VCO buffer output drive |

These are the TX-critical registers. The R820T2 has ~30 registers total;
only these 5 need to change from their receive-optimized defaults.

## When To Use

**Always** before running any Ether-VM mode that transmits:
- `--ofdm` / `--ofdm-file`
- `--vm` with ANTISYM, TICK, or TX instructions
- Any API program calling `qvm_ofdm_compute()`

Run once per SDR plug-in event. The registers are volatile — they reset on
power cycle or USB replug.

**Not needed** for receive-only modes: `--ether-gate`, `--ether-transfer`,
`--field`, `--ether-monitor`, `--bell` (loopback).

## Hardware Notes

- **RTL-SDR Blog V3/V4**: Some models have a hardware TX switch and may not
  need this boost. The boost is harmless if already TX-capable.
- **Generic R820T2 dongles**: These are the primary target. They ship as
  DVB-T receivers with no TX configuration at all.
- **0bda:2832 vs 0bda:2838**: Both use the R820T2. The device is probed
  at `lsusb -d 0bda:2838` but the code also detects 0bda:2832.

## Verification

After boosting, run a quick probe:

```sh
./sdr_ether 256 100e6 2048000 496 --vm
> INIT
> PROB
```

You should see non-zero bin powers across the spectrum, not just DC (bin 0).
If all power collapses to bin 0 within 2-3 TICK cycles, the boost may need
re-application or the register values need adjustment for your specific dongle
batches.
