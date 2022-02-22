
cd ihs
python ../../../../scripts/translate_coords.py inshore.ihs.50bins.norm.50kb.windows ../../ragtag/ragtag_output/ragtag.scaffolds.agp | awk '{OFS="\t";print $1,$2,$2+50000,$4,$5,$6,$7,$8}' > inshore.ihs.50bins.norm.50kb.windows.ragtag.txt
cd..

cd xpehh
for f in *.windows;do
	python ../../../../scripts/translate_coords.py ${f} ../../ragtag/ragtag_output/ragtag.scaffolds.agp | \
		awk '{OFS="\t";print $1,$2,$2+50000,$4,$5,$6,$7,$8,$9,$10}' > \
		${f}.ragtag.txt
done
cd ..

cd xpnsl
for f in *.windows;do
	python ../../../../scripts/translate_coords.py ${f} ../../ragtag/ragtag_output/ragtag.scaffolds.agp | \
		awk '{OFS="\t";print $1,$2,$2+50000,$4,$5,$6,$7,$8,$9,$10}' > \
		${f}.ragtag.txt
done
cd ..



python ../../../scripts/translate_coords.py pbs/plink2_noheader.pbs ../ragtag/ragtag_output/ragtag.scaffolds.agp >pbs/plink2.pbs_scaff.tsv


python ../../../scripts/translate_coords.py pbs/pbs_sweeps.tsv ../ragtag/ragtag_output/ragtag.scaffolds.agp > pbs/pbs_sweeps.scaff.tsv
