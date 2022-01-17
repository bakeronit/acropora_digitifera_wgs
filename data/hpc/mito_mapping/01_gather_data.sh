find /fast/shared/Acropora_digitifera_wgs_bamfile/ -name '*.bam' | xargs -I{} ln -s {} .
find /fast/shared/Acropora_digitifera_wgs_bamfile/ -name '*.bai' | xargs -I{} ln -s {} .
