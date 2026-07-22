#!/bin/bash
# boost.sh — Find RTL2838, set USB permissions, boost R820T2 TX registers

DEV=$(lsusb -d 0bda:2838 | head -1)
if [ -z "$DEV" ]; then
    echo "No RTL2838 found. Plug in the SDR."
    exit 1
fi

BUS=$(echo "$DEV" | awk '{print $2}' | sed 's/^0*//')
NUM=$(echo "$DEV" | awk '{print $4}' | tr -d ':' | sed 's/^0*//')
NUM3=$(printf '%03d' "$((10#$NUM))")
DEVPATH="/dev/bus/usb/00${BUS}/$NUM3"

echo "RTL2838 at Bus $BUS Device $NUM → $DEVPATH"

if [ ! -e "$DEVPATH" ]; then
    echo "Device node missing — mounting usbfs..."
    sudo mount -t usbfs none /proc/bus/usb 2>/dev/null
    DEVPATH="/proc/bus/usb/00${BUS}/$(printf '%03d' $NUM)"
fi

sudo chmod 666 "$DEVPATH" 2>/dev/null

if [ -x ./ether_boost ]; then
    echo "Running ether_boost..."
    ./ether_boost
else
    echo "Building ether_boost..."
    gcc -O2 ether_boost.c -o ether_boost && ./ether_boost
fi

echo ""
echo "R820T2 boosted. Verify: ls /dev/swradio*"
ls /dev/swradio* 2>/dev/null || echo "(no SDR device — run: sudo modprobe rtl2832_sdr)"
