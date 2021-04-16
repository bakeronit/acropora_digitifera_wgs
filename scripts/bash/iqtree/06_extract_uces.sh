# Requires bedtools v2.30.0
conda activate

bedtools getfasta -nameOnly -fi MI-1-4_merged_marked_regions.consensus.fa -bed aten_probes_merged.bed -fo MI-1-4.uces.fasta

bedtools getfasta -nameOnly -fi FI-1-3_merged_marked_regions.consensus.fa -bed aten_probes_merged.bed -fo FI-1-3.uces.fasta

bedtools getfasta -nameOnly -fi GCF_004143615.1_amil_sf_1.1_genomic.fna -bed amil_probes_merged.bed -fo amil.uces.fasta

bedtools getfasta -nameOnly -fi GCA_014634065.1_Adig_2.0_genomic.fna -bed adi_probes_merged.bed -fo adi.uces.fasta


for f in *aligned*.consensus.fa;do
	fn=${f%_aligned*}
	bedtools getfasta -nameOnly -fi $f -bed adi_probes_merged.bed -fo adi_${fn}.uces.fasta
done

