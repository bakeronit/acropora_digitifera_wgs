# This file is here simply to describe the process. 
# Actual script is /fast/shared/build_kraken_coral_symbionts2.sh on genomics2

#conda activate kraken1

DBNAME=kraken_coral_symbionts2


kraken-build --download-taxonomy --db $DBNAME

# Include standard bacteria in the database
#
kraken-build --download-library bacteria --db $DBNAME

find ./genome_krakenfa -name '*.fa' -print0 | \
    xargs -0 -I{} -n1 kraken-build --add-to-library {} --db $DBNAME

kraken-build --build --threads 16 --db $DBNAME

