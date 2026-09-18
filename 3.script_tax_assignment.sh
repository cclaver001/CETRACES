## For the taxonomic assignment we need:
# Samples_joined.fasta, Samples_joined.groups (the output files of the script 2.PrepareReads)

mothur

#create logfile with all manip that are going to be done
set.logfile(name=cetaceans_bob.log)

# Analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=Samples_joined.fasta, processors=12)
count.groups(group=Samples_joined.groups)

# Dereplicate
unique.seqs(fasta=Samples_joined.fasta)
summary.seqs(fasta=Samples_joined.unique.fasta, name=Samples_joined.names, processors=12)

# Align sequences against the reference database
alignmentlocation=/share/projects/OCEAN_eDNA/cetaceans_bob
align.seqs(fasta=Samples_joined.unique.fasta, reference=$alignmentlocation/Mitogenomes_cetacea_16S_newbarcoderegion.align, processors=12, flip=T)
summary.seqs(fasta=Samples_joined.unique.align, name=Samples_joined.names, processors=12)

# Remove sequences not covering the target region of 16S and with ambiguous bases
screen.seqs(fasta=Samples_joined.unique.align, name=Samples_joined.names, group=Samples_joined.groups, minlength=140, maxlength=230, maxambig=0, processors=12)
summary.seqs(fasta=Samples_joined.unique.good.align, name=Samples_joined.good.names, processors=12)
count.groups(group=Samples_joined.good.groups)

# Remove columns that contain gap characters
filter.seqs(fasta=Samples_joined.unique.good.align, vertical=T, processors=12) #elimina columnas del alineamiento donde todas las secuencias tienen gaps
unique.seqs(fasta=Samples_joined.unique.good.filter.fasta, name=Samples_joined.good.names)
summary.seqs(fasta=Samples_joined.unique.good.filter.unique.fasta,name=Samples_joined.unique.good.filter.names,processors=12)

# Remove chimeras
chimera.uchime(fasta=Samples_joined.unique.good.filter.unique.fasta, name=Samples_joined.unique.good.filter.names, group=Samples_joined.good.groups, processors=12)
remove.seqs(accnos=Samples_joined.unique.good.filter.unique.denovo.uchime.accnos, fasta=Samples_joined.unique.good.filter.unique.fasta, name=Samples_joined.unique.good.filter.names, group=Samples_joined.good.groups, dups=T)
summary.seqs(fasta=Samples_joined.unique.good.filter.unique.pick.fasta, name=Samples_joined.unique.good.filter.pick.names, processors=12)
count.groups(group=Samples_joined.good.pick.groups)

# Rename files
system(cp Samples_joined.unique.good.filter.unique.pick.fasta Samples_joined_all.fasta)
system(cp Samples_joined.unique.good.filter.pick.names Samples_joined_all.names)
system(cp Samples_joined.good.pick.groups Samples_joined_all.groups)

# Create count_table with the new files
make.table(name=Samples_joined_all.names, group=Samples_joined_all.groups)

quit()

# Clustering by SWARM
Rscript ../GIT/3.script_prep_SWARM.R Samples_joined_all.names Samples_joined_all.fasta # modify files for swarm
swarm Samples_joined_all_swarm.fasta -d 1 -t 12 -f -l Samples_joined_all_swarm_d1.log -r -o Samples_joined_all_swarm_d1.out -s Samples_joined_all_swarm_d1.stats -w Samples_joined_all_swarm_d1_rep.fasta
sed -E  s/\(_[0-9]*,\)/","/g Samples_joined_all_swarm_d1.out | sed -E s/\(_[0-9]*\\t\)/"\\t"/g | sed -E  s/\(_[0-9]*$\)//g | sed s/swarm/swarm_1/g > Samples_joined_all_swarm_d1.list # make list file
mothur "#make.shared(list=Samples_joined_all_swarm_d1.list, count=Samples_joined_all.count_table, label=swarm_1)" # make shared file
mothur "#remove.rare(shared=Samples_joined_all_swarm_d1.shared, nseqs=1)" # remove singletons
mv Samples_joined_all_swarm_d1.swarm_1.pick.shared Samples_joined_all_swarm_d1.noSingletons.shared # rename
awk '/^>/{printf ">ASV%03d\n",++i; next}{print}' Samples_joined_all_swarm_d1_rep.fasta > Samples_joined_all_swarm_d1_rep_otuName.fasta #cuidado con los ceros!!!!!!
# Remove singletons in fasta file also, based on singletons removed in the shared file in the previous step
cat  Samples_joined_all_swarm_d1.noSingletons.shared | awk 'NR==1,NR==1' | cut -f 1,2,3 --complement | sed 's/\t/\n/g' >  Samples_joined_all_swarm_d1.noSingletons.OTUs_temp
awk 'NR==FNR { keep[$1]; next } /^>/ { id = substr($1,2) } id in keep { print; getline; print; delete keep[id] }'  Samples_joined_all_swarm_d1.noSingletons.OTUs_temp Samples_joined_all_swarm_d1_rep_otuName.fasta > Samples_joined_all_swarm_d1_rep_otuName.noSingletons.fasta
rm -f Samples_joined_all_swarm_d1.noSingletons.OTUs_temp

## LULU (post-clustering correction)
makeblastdb -in Samples_joined_all_swarm_d1_rep_otuName.noSingletons.fasta -parse_seqids -dbtype nucl
blastn -db Samples_joined_all_swarm_d1_rep_otuName.noSingletons.fasta -outfmt '6 qseqid sseqid pident' -out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query Samples_joined_all_swarm_d1_rep_otuName.noSingletons.fasta -num_threads 16
cut -f 1,3 --complement Samples_joined_all_swarm_d1.noSingletons.shared > Samples_joined_all_swarm_d1.noSingletons.table.tsv
cp Samples_joined_all_swarm_d1.noSingletons.table.tsv SWARM_table_clean.tsv
Rscript ../GIT/3.script_prep_LULU.R # output SWARM_table_curated.tsv
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < SWARM_table_curated.tsv > SWARM_table_curated_t.tsv #traspose file
# Remove SWARMS removed by LULU in fasta
cut -f 1 SWARM_table_curated_t.tsv | tail -n +2 > Samples_joined_all_swarm_d1.noSingletons.curated.OTUs_temp
awk 'NR==FNR { keep[$1]; next } /^>/ { id = substr($1,2) } id in keep { print; getline; print; delete keep[id] }' Samples_joined_all_swarm_d1.noSingletons.curated.OTUs_temp Samples_joined_all_swarm_d1_rep_otuName.noSingletons.fasta > Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.fasta

# Remove gaps from sequs to perform the tax assignment
mothur "#degap.seqs(fasta= Samples_joined_all.fasta)" #without SWARM and LULU
mothur "#degap.seqs(fasta= Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.fasta)" # output Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.ng.fasta

#TAXONOMIC ASSIGNMENT
DBlocation=/share/projects/OCEAN_eDNA/cetaceans_bob
DBprefix=Cetacea_16S_REFDB_newbarcoderegion
mothur "#classify.seqs(fasta=Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.ng.fasta, template=$DBlocation/$DBprefix.fasta, taxonomy=$DBlocation/$DBprefix.tax, name=Samples_joined_all.names, group=Samples_joined_all.groups, method=wang, cutoff=80, processors=12)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.ng.Cetacea_16S_REFDB_newbarcoderegion.wang.taxonomy | sort -k1 | sed '1s/^/OTU\tTaxonomy\n/' > Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.ng.Cetacea_16S_REFDB_newbarcoderegion.wang_corr.taxonomy
paste SWARM_table_curated_t.tsv Samples_joined_all_swarm_d1_rep_otuName.noSingletons.curated.ng.Cetacea_16S_REFDB_newbarcoderegion.wang_corr.taxonomy > OTUtable_CetaceaJUVENA.txt
# to merge samples in phylotypes (instead of SWARM)
mothur "#phylotype(taxonomy=Samples_joined_all.ng.$DBprefix.wang.taxonomy)"
mothur "#make.shared(list=Samples_joined_all.ng.$DBprefix.wang.tx.list, count=Samples_joined_all.count_table, label=1)"
mothur "#classify.otu(list=Samples_joined_all.ng.$DBprefix.wang.tx.list, count=Samples_joined_all.count_table, taxonomy=Samples_joined_all.ng.$DBprefix.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g Samples_joined_all.ng.$DBprefix.wang.tx.1.cons.taxonomy > Samples_joined_all.ng.$DBprefix.wang.tx.1.cons_corr.taxonomy
#merge information and create output files
awk '{$1=$3=""; print $0}' Samples_joined_all.ng.$DBprefix.wang.tx.shared > file1 #remove label and phylo_count columns
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < file1 > file2 #transpose file
paste file2 Samples_joined_all.ng.$DBprefix.wang.tx.1.cons_corr.taxonomy > cetaceans_bob_ng_Phylotable_70.txt

# BLAST ALIGNMENT
# 1. Create the BLAST local reference databse. # First, we have included the names of the species to which each sequence corresponds to in a new fasta file. Example: >NC_001601.1_Balaenoptera_musculus AAATCACAACCTTAAACCACCAAGG
#/share/projects/OCEAN_eDNA/cetaceans_bob/tax_assign_BLAST/ just to remember the folder
makeblastdb -in /share/projects/OCEAN_eDNA/cetaceans_bob/tax_assign_BLAST/Cetacea_16S_REFDB_newbarcoderegion_withnames.fasta  -dbtype nucl  -out Cetacea_16S_with_names

# 2. Run BLAST
# -max_target_seqs número de hits que se reportan
# dentro del outfmt: qseqid (ID del read), sseqid (ID de la secuencia de referencia con la que matchea), pident (% de identidad del alineamiento), length (longitud del alineamiento), mismatch (número de mismatches), qstart/qend (posición de inicio/final en la consulta), sstart/end (inicio/final en la secuencia sujeto)
#-max_hsps: devolver hasta X sub‑alineamientos contra la MISMA secuencia de referencia
blastn -query ../prepare_reads/Samples_joined_all.fasta -db Cetacea_16S_with_names -max_hsps 1 -max_target_seqs 3 -evalue 0.0001 -perc_identity 95 -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames scomnames" -out Samples_joined_all.BLAST_with_tax_3_hits.txt
blastn -query ../prepare_reads/Samples_joined_all.fasta -db Cetacea_16S_with_names -max_hsps 1 -max_target_seqs 1 -evalue 0.0001 -perc_identity 95 -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames scomnames" -out Samples_joined_all.BLAST_with_tax_1_hit.txt

paste ../prepare_reads/Samples_joined_all.count_table Samples_joined_all.BLAST_with_tax_1_hit.txt > BLAST-results-table.txt


awk -F'\t' '{print $1"\t"$2";"}' Samples_joined_all.BLAST_with_tax_1_hit.txt | sed '1s/^/OTU\tTaxonomy\n/' > BLAST.taxonomy
# merge taxonomy with count_table

tail -n +2 ../prepare_reads/Samples_joined_all.count_table | sort -t $'\t' -k1,1 > file1_sorted.tsv
sort -t $'\t' -k1,1 Samples_joined_all.count_table_with_taxonomy_BLAST.txt > file2_sorted.tsv
join -t $'\t' -a1 -a2 -e "NA" -o auto file1_sorted.tsv file2_sorted.tsv > merged_body.tsv
