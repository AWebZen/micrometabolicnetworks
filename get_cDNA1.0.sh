if [[ $# -lt 2 ]] ; then
    echo 'Missing arguments. Needs an assembly file first and organism name next.'
    exit 1
fi

org=`echo "$2"|sed 's/_/\ /g'`
org_modif=`echo "$org"|sed 's/\ /\.\*/g'`

firsthit=`grep -m 1 "$org_modif" $1`
if [[ -z $firsthit ]] #No hits in ensembl file
then
	firsthit=`grep -m 1 "${org##* }" $1`
	#echo -e "NOT FOUND ON FIRST TRY, CHECK IF OK (if empty, no hit on second try):\n$firsthit"
fi
if [[ -z $firsthit ]] #No hits in ensembl file in second try
then
	echo "No hit for species $org"
    address="None"
else
	collection_core=`echo "$firsthit"|awk -F"\t" '{print $13}'` #collection for url address
	collection=${collection_core%%_core*} #cropped
	species=`echo "$firsthit"|awk -F"\t" '{print $2}'` #correct species names for url address
	address="ftp://ftp.ensemblgenomes.org/pub/release-38/bacteria//fasta/$collection/$species/cdna/*cdna.all.fa.gz"
fi

wget $address -P genomes_cdna
mv genomes_cdna/${address##*/} genomes_cdna/$2"_cds_from_genomic.fna.gz"
gzip -d genomes_cdna/*_cds_from_genomic.fna.gz
rm -f genomes_cdna/*_cds_from_genomic.fna.gz


