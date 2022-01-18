# This script is for updating shared data in the cloudstor repository
# It is designed for collaborators with direct access to the cloudstor files
# The idea is that this will allow fast updates to the data because it only updates modified files
#

# The files to be synced are only those listed in data.list and data_large.list
# 
# data.list > Include any data files not checked into git (typically anything other than code) but <50Mb
# data_large.list > Is for very large files that we expect few users will actually want to download (>50Mb)
#

#
# The rclone parts of this require some setup before they will work
# See https://support.aarnet.edu.au/hc/en-us/articles/115007168507 for details
# 
#


rclone copy -cu -P --files-from jia-data.list . CloudStor:/acropora_digitifera_wgs/data_github
#rclone copy -cu -P --files-from data_large.list . CloudStor:/acropora_digitifera_wgs/data_github
