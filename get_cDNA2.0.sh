if [[ $# -lt 3 ]] ; then
    echo 'Missing arguments. Needs an assembly file first, organism name and organism taxonomy code next.'
    exit 1
fi

address=`grep -m 1 "	$3	" $1|awk -F"\t" '{print $20}'|sed 's/\(GCA_.*\)$/\1\/\1_cds_from_genomic.fna.gz/g'`
wget $address -P genomes_cdna
mv genomes_cdna/${address##*/} genomes_cdna/$2"_cds_from_genomic.fna.gz"
gzip -d genomes_cdna/*_cds_from_genomic.fna.gz
rm -f genomes_cdna/*_cds_from_genomic.fna.gz


