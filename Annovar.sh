#!/bin/bash

##ANNOVAR
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/

perl annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/

perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/

perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct humandb/

perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/

perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all humandb/

#Create the input files from the vcf files
perl convert2annovar.pl -format vcf4 -allsample -withfreq /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.biallelic.vcf  > /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.biallelic.avinput

#Create the annotation
perl table_annovar.pl /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.biallelic.avinput humandb/ -buildver hg19 -out myanno.joint.output.snp.indel.recal.avinput -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -csvout