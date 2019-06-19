outf=species_classified_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done


outf=species_classified_genbank2018.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_genbank2018/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_genbank2019.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_genbank2019/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_2018_unmapped_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_2018_unmapped_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_2019_unmapped_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_2019_unmapped_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done


## META1k
outf=species_classified_genbank2018.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_genbank_hq/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_genbank2019.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_jan2019/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_2018_unmapped_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_2018_unmapped_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

outf=species_classified_2019_unmapped_segata.txt
printf "%s\t%s\t%s\t%s\n" "file" "species" "total" "frac" > $outf
for i in kraken2_classification_2019_unmapped_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\n" $(basename $i) $species $total $frac >> $outf
done

# fionas south africa samples
outf=species_classified_feb2019.txt
printf "%s\t%s\t%s\t%s\t%s\n" "file" "species" "root" "total" "frac" > $outf
for i in kraken2_classification_feb2019/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    root=$(awk '$4=="R" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\t%s\n" $(basename $i) $species $root $total $frac >> $outf
done

outf=species_classified_feb2019_bracken.txt
printf "%s\t%s\t%s\t%s\t%s\n" "file" "species" "root" "total" "frac" > $outf
for i in kraken2_classification_feb2019/classification/*.krak_bracken.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    root=$(awk '$4=="R" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\t%s\n" $(basename $i) $species $root $total $frac >> $outf
done

outf=species_classified_feb2019_unmapped_segata.txt
printf "%s\t%s\t%s\t%s\t%s\n" "file" "species" "root" "total" "frac" > $outf
for i in kraken2_classification_feb2019_unmapped_segata/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    root=$(awk '$4=="R" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\t%s\n" $(basename $i) $species $root $total $frac >> $outf
done

outf=species_classified_feb2019_unmapped_nayfach.txt
printf "%s\t%s\t%s\t%s\t%s\n" "file" "species" "root" "total" "frac" > $outf
for i in kraken2_classification_feb2019_unmapped_nayfach/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    root=$(awk '$4=="R" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\t%s\n" $(basename $i) $species $root $total $frac >> $outf
done

outf=species_classified_feb2019_unmapped_almeida.txt
printf "%s\t%s\t%s\t%s\t%s\n" "file" "species" "root" "total" "frac" > $outf
for i in kraken2_classification_feb2019_unmapped_almeida/classification/*.krak.report; do
    species=$(awk '$4=="S" { s+=$2 } END { print s }' $i) 
    root=$(awk '$4=="R" { s+=$2 } END { print s }' $i) 
    total=$(head -n2 $i | awk '{ s+=$2 } END { print s }')
    frac=$( echo "scale=4; $species / $total * 100" | bc)
    printf "%s\t%s\t%s\t%s\t%s\n" $(basename $i) $species $root $total $frac >> $outf
done
