#!/bin/bash

echo "Running ntRoot demo..."
echo "Please ensure that ntRoot is on your PATH"

set -eux -o pipefail

echo "Running ntRoot reads demo..."
if [ ! -f ERR3239334.chr21_1.fq.gz ]; then
	wget https://www.bcgsc.ca/downloads/btl/ntroot/reads_demo/ERR3239334.chr21_1.fq.gz
fi
if [ ! -f ERR3239334.chr21_2.fq.gz ]; then
	wget https://www.bcgsc.ca/downloads/btl/ntroot/reads_demo/ERR3239334.chr21_2.fq.gz
fi

ntroot --reference chr21.fa --reads ERR3239334.chr21_ -k 55 -l pop-spec-snp_chr21.vcf.gz --lai

prediction=$(cat ERR3239334.chr21__ntedit_k55_variants.vcf_ancestry-predictions_tile5000000.tsv | grep -v "#" | awk '{print $1}' |head -n 2 |tail -n 1)

if [ ${prediction} == "EUR" ]; then
	echo "ntRoot reads test successful!"
else 
	echo "ntRoot reads test failed.. Please check your installation."
	exit 1
fi

echo "Running ntRoot genome demo..."

ntroot --reference chr21.fa --genome HuRef.chr21.fa -k 55 -l pop-spec-snp_chr21.vcf.gz

prediction=$(cat HuRef.chr21_ntedit_k55_variants.vcf_ancestry-predictions_tile5000000.tsv | grep -v "#" | awk '{print $1}' |head -n 2 |tail -n 1)
if [ ${prediction} == "EUR" ]; then
	echo "ntRoot genome test successful!"
else 
	echo "ntRoot genome test failed.. Please check your installation."
	exit 1
fi

echo "Running ntRoot input VCF demo..."

ntroot --reference chr21.fa --custom_vcf HG002.chr21.vcf.gz -l pop-spec-snp_chr21.vcf.gz
prediction=$(cat HG002.chr21.vcf.gz.cross-ref.vcf_ancestry-predictions_tile5000000.tsv | grep -v "#" | awk '{print $1}' |head -n 2 |tail -n 1)
if [ ${prediction} == "EUR" ]; then
	echo "ntRoot input VCF test successful!"
else 
	echo "ntRoot input VCF test failed.. Please check your installation."
	exit 1
fi

echo "Running ntRoot --exome reads demo with --exome_bed..."

if [ ! -f HG00864_ERR050736.chr20-21.fastq.gz ]; then
	wget https://www.bcgsc.ca/downloads/btl/ntroot/reads_demo/HG00864_ERR050736.chr20-21.fastq.gz
fi
if [ ! -f HG00864_ERR050737.chr20-21.fastq.gz ]; then
	wget https://www.bcgsc.ca/downloads/btl/ntroot/reads_demo/HG00864_ERR050737.chr20-21.fastq.gz
fi

ntroot --reference chr20-21.fa.gz --reads HG00864 -k 55 -l pop-spec-snp_chr20-21.vcf.gz --exome --exome_bed exome_targets.bed
prediction=$(cat HG00864_ntedit_k55_exome_variants.vcf_ancestry-predictions_tile5000000.tsv | grep -v "#" | awk '{print $1}' |head -n 2 |tail -n 1)
if [ ${prediction} == "EAS" ]; then
	echo "ntRoot --exome reads test successful!"
else 
	echo "ntRoot --exome reads test failed.. Please check your installation."
	exit 1
fi

echo "Running ntRoot --exome demo with --masked (input reference already masked based on exon target coordinates)"
ntroot --exome --reference masked_chr20-21.fa.gz --reads HG00864 -k 55 -l pop-spec-snp_chr20-21.vcf.gz --masked --force
prediction=$(cat HG00864_ntedit_k55_exome_variants.vcf_ancestry-predictions_tile5000000.tsv | grep -v "#" | awk '{print $1}' |head -n 2 |tail -n 1)
if [ ${prediction} == "EAS" ]; then
	echo "ntRoot --exome masked test successful!"
else 
	echo "ntRoot --exome masked test failed.. Please check your installation."
	exit 1
fi


echo "Done ntRoot tests!"
