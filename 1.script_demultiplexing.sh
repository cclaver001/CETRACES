### PIPELINE FOR DEMULTIPLEXING -Tagsteady protocol- ###
mix=(`cut -f 1 /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/Mix_path_to_R1_R2_files.txt`)
R1=(`cut -f 2 /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/Mix_path_to_R1_R2_files.txt`)
R2=(`cut -f 3 /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/Mix_path_to_R1_R2_files.txt`)
x=(`wc -l /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/Mix_path_to_R1_R2_files.txt`)
y=$((x-1))

echo "Mix,MixRaw,Sample,RawReads_R1_step1,RawReads_R1_step2,RawReads_R2_step1,RawReads_R2_step2" > reads_stats_samples.txt

for i in `seq 0 $y`; do
  gzip -dc < ${R1[i]} > ${mix[i]}_R1.fastq
  gzip -dc < ${R2[i]} > ${mix[i]}_R2.fastq
  r="$(grep -c '^@M0' ${mix[i]}_R1.fastq)"
  echo "" > ${mix[i]}.log

  sample=(`cut -f 1 /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/${mix[i]}_barcodes.txt`)
  Tags=(`cut -f 2 /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/${mix[i]}_barcodes.txt`)
  a=(`wc -l /share/projects/OCEAN_eDNA/cetaceans_bob/GIT/${mix[i]}_barcodes.txt`)
  b=$((a-1))

  for j in `seq 0 $b`; do
    cutadapt -g ${Tags[j]} -G ${Tags[j]} --action=none --discard-untrimmed -e 0 -O 6 -o ${sample[j]}_${mix[i]}_R1_temp.fastq -p ${sample[j]}_${mix[i]}_R2_temp.fastq ${mix[i]}_R1.fastq ${mix[i]}_R2.fastq >> ${mix[i]}.log
    cutadapt -g ${Tags[j]} -G ${Tags[j]} --discard-untrimmed -e 0 -O 6 -o ${sample[j]}_${mix[i]}_R2.fastq -p ${sample[j]}_${mix[i]}_R1.fastq ${sample[j]}_${mix[i]}_R2_temp.fastq ${sample[j]}_${mix[i]}_R1_temp.fastq >> ${mix[i]}.log
# -e es cero porque no queremos ninguna variación
# action none es para que no elimine la secuencia identificada
#los fastas temporales son para los que sólo se ha identificado uno de los dos primers y necesitan una segunda pasada
    echo -n ${mix[i]} >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    echo -n $r >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    echo -n ${sample[j]} >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    s="$(grep -c '^@M0' ${sample[j]}_${mix[i]}_R1_temp.fastq)"
    echo -n $s >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    t="$(grep -c '^@M0' ${sample[j]}_${mix[i]}_R1.fastq)"
    echo -n $t >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    u="$(grep -c '^@M0' ${sample[j]}_${mix[i]}_R2_temp.fastq)"
    echo -n $u >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    v="$(grep -c '^@M0' ${sample[j]}_${mix[i]}_R2.fastq)"
    echo -n $v >> reads_stats_samples.txt
    echo -e "" >> reads_stats_samples.txt
  done
rm -f *_temp.fastq
done
