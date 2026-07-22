/*
 * ether_boost — User-space R820T2 register poke via USB control transfer
 *
 * One-time: sudo chmod 666 /dev/bus/usb/001/016
 * Then:     ./ether_boost
 *
 * Build: gcc -O2 ether_boost.c -o ether_boost
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <linux/usbdevice_fs.h>

#define USB_TYPE_VENDOR (0x02 << 5)
#define USB_DIR_OUT 0

static int r820t2_write(int fd, uint8_t reg, uint8_t val) {
    struct usbdevfs_ctrltransfer ctrl;
    uint8_t buf[3];
    memset(&ctrl, 0, sizeof(ctrl));
    ctrl.bRequestType = USB_TYPE_VENDOR | USB_DIR_OUT;
    ctrl.bRequest = 0x02;
    ctrl.wValue = 0;
    ctrl.wIndex = 0;
    ctrl.wLength = 3;
    ctrl.data = buf;
    ctrl.timeout = 500;
    buf[0] = 0x34; /* R820T2 I2C addr << 1 */
    buf[1] = reg;
    buf[2] = val;
    return ioctl(fd, USBDEVFS_CONTROL, &ctrl);
}

int main(int argc, char **argv) {
    int fd = open("/dev/bus/usb/001/016", O_RDWR);
    if (fd < 0) {
        /* Try to find the device */
        for (int n = 1; n < 32; n++) {
            char path[64];
            snprintf(path, 64, "/dev/bus/usb/001/%03d", n);
            fd = open(path, O_RDWR);
            if (fd >= 0) break;
        }
    }
    if (fd < 0) {
        fprintf(stderr, "Cannot open USB device.\n");
        fprintf(stderr, "Run: lsusb -d 0bda:2838\n");
        fprintf(stderr, "Then: sudo chmod 666 /dev/bus/usb/XXX/YYY\n");
        return 1;
    }

    printf("Boosting R820T2 registers...\n");
    r820t2_write(fd, 0x05, 0x00);  /* LNA bypass */
    r820t2_write(fd, 0x06, 0x01);  /* LNA bypass enable */
    r820t2_write(fd, 0x07, 0x1F);  /* max mixer bias */
    r820t2_write(fd, 0x0A, 0x60);  /* max VCO amplitude */
    r820t2_write(fd, 0x10, 0x0F);  /* VCO buffer amp */
    printf("Done — R820T2 boosted (LNA bypass, max VCO/mixer).\n");
    close(fd);
    return 0;
}
