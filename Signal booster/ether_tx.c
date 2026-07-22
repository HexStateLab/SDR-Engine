/*
 * ether_tx.c — Kernel module: boost RTL-SDR TX via R820T2 I2C
 *
 * Registers as a USB driver for RTL2832/2838.
 * On probe, creates /sys/kernel/ether_tx/gain (0-100).
 * Writing a value sends I2C commands through the RTL2832U.
 *
 * Build: make
 * Usage:
 *   sudo insmod ether_tx.ko
 *   echo 100 > /sys/kernel/ether_tx/gain
 *   sudo rmmod ether_tx
 */

#include <linux/module.h>
#include <linux/kernel.h>
#include <linux/init.h>
#include <linux/usb.h>
#include <linux/kobject.h>
#include <linux/sysfs.h>
#include <linux/slab.h>

MODULE_LICENSE("GPL");
MODULE_AUTHOR("SDRQudit");
MODULE_DESCRIPTION("RTL-SDR TX power boost via R820T2 registers");

static struct usb_device *rtl_dev;
static int tx_gain = 50;
static struct kobject *ether_kobj;
static struct kobj_attribute gain_attr;

#define RTL2832U_REQ_I2C 0x02

static int r820t2_write(struct usb_device *dev, u8 reg, u8 val)
{
    u8 buf[3];
    buf[0] = 0x34; /* R820T2 I2C addr << 1 */
    buf[1] = reg;
    buf[2] = val;

    return usb_control_msg(dev, usb_sndctrlpipe(dev,0),
        RTL2832U_REQ_I2C, USB_TYPE_VENDOR | USB_DIR_OUT,
        0, 0, buf, 3, 1000);
}

static void r820t2_boost(struct usb_device *dev, int level)
{
    u8 lna = (level > 80) ? 0x00 : 0x20;
    u8 mixer = 0x0F;
    u8 vco = 0x30;
    r820t2_write(dev, 0x05, lna);
    r820t2_write(dev, 0x07, mixer);
    r820t2_write(dev, 0x0A, vco);
    pr_info("ether_tx: LNA=%02x MIX=%02x VCO=%02x level=%d\n",lna,mixer,vco,level);
}

static ssize_t gain_show(struct kobject *k, struct kobj_attribute *a, char *b){return sprintf(b,"%d\n",tx_gain);}
static ssize_t gain_store(struct kobject *k, struct kobj_attribute *a, const char *b, size_t n){
    int v;if(kstrtoint(b,10,&v)||v<0||v>100)return -EINVAL;
    tx_gain=v;r820t2_boost(rtl_dev,tx_gain);return n;}

static int ether_probe(struct usb_interface *intf, const struct usb_device_id *id)
{
    rtl_dev = interface_to_usbdev(intf);
    usb_get_dev(rtl_dev);
    pr_info("ether_tx: probed %04x:%04x\n",
        le16_to_cpu(rtl_dev->descriptor.idVendor),
        le16_to_cpu(rtl_dev->descriptor.idProduct));

    ether_kobj = kobject_create_and_add("ether_tx", kernel_kobj);
    if(!ether_kobj) return -ENOMEM;
    gain_attr = (struct kobj_attribute)__ATTR(gain,0664,gain_show,gain_store);
    sysfs_create_file(ether_kobj,&gain_attr.attr);

    r820t2_boost(rtl_dev, tx_gain);
    return 0;
}

static void ether_disconnect(struct usb_interface *intf)
{
    if(rtl_dev){r820t2_write(rtl_dev,0x05,0x20);usb_put_dev(rtl_dev);rtl_dev=NULL;}
    if(ether_kobj){sysfs_remove_file(ether_kobj,&gain_attr.attr);kobject_put(ether_kobj);}
    pr_info("ether_tx: disconnected\n");
}

static struct usb_device_id ether_table[] = {
    {USB_DEVICE(0x0bda,0x2832)},
    {USB_DEVICE(0x0bda,0x2838)},
    {}
};
MODULE_DEVICE_TABLE(usb, ether_table);

static struct usb_driver ether_driver = {
    .name = "ether_tx",
    .id_table = ether_table,
    .probe = ether_probe,
    .disconnect = ether_disconnect,
};

module_usb_driver(ether_driver);
