# This script is for building user friendly data packages for end-users
# By packaging all files into a small and large tarball we allow users without 
# direct cloudstor access to download the data via a public link which we publish
# in the README

# For this to work it is important that the file listings data.list and data_large.list
# are maintained.  They should contain all files referenced by the RMarkdown. 
#
# data.list > Include any data files not checked into git (typically anything other than code) but <50Mb
# data_large.list > Is for very large files that we expect few users will actually want to download (>50Mb)
#

# The rclone parts of this require some setup before they will work
# See https://support.aarnet.edu.au/hc/en-us/articles/115007168507 for details
# 
# You will also need to ensure that the path Shared/wa_digitifera exists
#

DIFF=$(diff <(cat data.list | grep -v '^#' | grep -v '^$') <(tar -tf data.tgz))
DIFFLARGE=$(diff <(cat data_large.list | grep -v '^#' | grep -v '^$') <(tar -tf data_large.tgz))
if [ "$DIFF" != "" ]; then
	echo "Rebuilding data.tgz"
	tar -zcvf data.tgz -T <(cat data.list | grep -v '^#')
	# Upload to cloudstor if needed
	rclone copy --progress --no-traverse data.tgz CloudStor:/Shared/wa_digitifera/
	echo "Done uploading data.tgz"
fi


if ! [ -e data_large.tgz ] || [ "$DIFFLARGE" != "" ]; then
	echo "Rebuilding data_large.tgz"
	tar -zcvf data_large.tgz -T <(cat data_large.list | grep -v '^#')
	# Upload to cloudstor if needed
	rclone copy --progress --no-traverse data_large.tgz CloudStor:/Shared/wa_digitifera/
	echo "Done uploading data_large.tgz"
fi

