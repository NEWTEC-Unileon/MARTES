
#!/bin/bash
source /home/josecobo/anaconda3/etc/profile.d/conda.sh
conda activate platon
echo "----> Runing platon for plasmid detection"
platon --db /datos/DATABASE/platon_db/db/ --output ARGcontigs_platon --threads 32 ARG_contigs.fna
plasmidfinder.py -i ARG_contigs.fna -o plasmidfinder -x
conda deactivate

kraken2 --db /datos/DATABASE/PlusPF_20251015/ --use-names --threads 32 --output ARGcontigs_kraken2_PlusPF_20251015.txt ARG_contigs.fna
kraken2 --db /datos/DATABASE/PlusPF_20251015/ --use-names --threads 32 --confidence 0.1 --output ARGcontigs_kraken2_PlusPF_20251015_01.txt --report ARGcontigs_kraken2_report_01.txt ARG_contigs.fna
grep -P "\tG\t" ARGcontigs_kraken2*_report*.txt | cut -f 6 | sed "s/ //g" > genera_names.txt

conda activate genomad
genomad end-to-end ARG_contigs.fna ARGcontigs_genomad /datos/DATABASE/genomad_db/
conda deactivate

####		nano ARGcontigs_genomad/ARG_contigs_summary/*summary.tsv

# no funciona en delisoil (problemas con pandas, numpy y su puta madre (lo corro en master)
integron_finder --cpu 1 --outdir ARGcontigs_integronfinder ARG_contigs.fna

conda activate MobileElementFinder
mefinder find -c ARG_contigs.fna -t 32 mefinder
conda deactivate
