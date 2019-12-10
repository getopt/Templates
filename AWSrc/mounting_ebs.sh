

# sudo mkfs -t ext4 /dev/xvdf # reformat if needed, e.g. a fresh volume

sudo mkdir /data/
sudo mount /dev/xvdf /data/

sudo chown ec2-user:ec2-user data/




#
# /etc/fstab approach
#

# After adding to /etc/fstab line:
# /dev/xvdba   /data       ext4    defaults,nofail 0   0
# 

sudo mount -a
sudo chown ec2-user:ec2-user data/
