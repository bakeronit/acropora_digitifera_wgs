# This script is for updating shared data in the cloudstor repository
# It is designed for collaborators with direct access to the cloudstor files
#
#
# See https://support.aarnet.edu.au/hc/en-us/articles/115007168507 for details of setting up rclone
# 
#
rclone copy -P -u CloudStor:/Shared/data_github/data/hpc data/hpc/

rclone copy -P -u CloudStor:/Shared/data_github/data/maps data/maps/
