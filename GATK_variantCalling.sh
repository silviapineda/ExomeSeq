####This scripts is run in the server where the BAM files are

####GATK to call GVCF files
FILES=/home/pinedasans/ExomeSeq/BAM/*.bam


##This is to create the index (.bai) for each bam file
for f in $FILES
do
    echo $f
    samtools index $f
    var=$(basename $f)
    echo $var
    
    ##Running the HaplotypeCaller to generate the gvcf files
    java -jar GenomeAnalysisTK.jar -nct 20 -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -T HaplotypeCaller -I $f --emitRefConfidenc
    e GVCF --dbsnp /home/pinedasans/ExomeSeq/BAM/dbsnp.vcf -o /home/pinedasans/ExomeSeq/VCF/$var.raw.sn
    ps.indels.g.vcf
done

for f in $FILES
do
    var=$(basename $f)
    variantParameters="$variantParameters --variant /home/pinedasans/ExomeSeq/VCF/$var.raw.snps.indels.g.vcf"
done

echo start

#Obtain the genotypes and merge the samples
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa $variantParameters -o /home/pinedas
ans/ExomeSeq/VCF/joint.output.vcf

echo EndGenotypes

##To delete genotypes with minDP 8 and minGQ 20
##To run vcftools first: export PATH=${PATH}:/home/pinedasans/programs/vcftools_0.1.13/bin/
#When I applied this filer the variant recalibrator is not working.
#vcftools --minDP 8 --minGQ 20 --vcf /home/pinedasans/ExomeSeq/VCF/joint.output.vcf --recode --out /home/pinedasans/ExomeSeq/VCF/joint.output.DP.GQ.filtered.vcf


####Recalibration Best practice GATK
#https://software.broadinstitute.org/gatk/guide/article?id=2805
#https://software.broadinstitute.org/gatk/guide/article?id=39

##Recalibrate SNPs
java -jar GenomeAnalysisTK.jar -T VariantRecalibrator -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -input /home/pinedasans/ExomeSeq/VCF/joint.output.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/pinedasans/ExomeSeq/VCF/gVCF/hapmap_3.3.b37.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 /home/pinedasans/ExomeSeq/VCF/gVCF/1000G_omni2.5.b37.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/pinedasans/ExomeSeq/VCF/gVCF/1000G_phase1.snps.high_confidence.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/pinedasans/ExomeSeq/VCF/gVCF/dbsnp_138.b37.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP -tranche 99.9 -recalFile /home/pinedasans/ExomeSeq/VCF/output.recal -tranchesFile /home/pinedasans/ExomeSeq/VCF/output.tranches -rscriptFile /home/pinedasans/ExomeSeq/VCF/output.plots.R

echo EndRecalibrateSNP

#Apply recalibration from SNPs
java -jar GenomeAnalysisTK.jar -T ApplyRecalibration -input /home/pinedasans/ExomeSeq/VCF/joint.output.vcf -o /home/pinedasans/ExomeSeq/VCF/joint.output.snp.recal.vcf  -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa --ts_filter_level 99.9 -tranchesFile /home/pinedasans/ExomeSeq/VCF/output.tranches -recalFile /home/pinedasans/ExomeSeq/VCF/output.recal -mode SNP

echo EndApplyRecalibrateSNP

###Recalibrate Indels
java -jar GenomeAnalysisTK.jar -T VariantRecalibrator -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -input /home/pinedasans/ExomeSeq/VCF/joint.output.snp.recal.vcf -resource:mills,known=false,training=true,truth=true,prior=12.0 /home/pinedasans/ExomeSeq/VCF/gVCF/Mills_and_1000G_gold_standard.indels.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/pinedasans/ExomeSeq/VCF/gVCF/dbsnp_138.b37.vcf -an MQRankSum -an ReadPosRankSum -an FS -mode INDEL -tranche 99.9 --maxGaussians 4 -recalFile /home/pinedasans/ExomeSeq/VCF/output.indel.recal -tranchesFile /home/pinedasans/ExomeSeq/VCF/output.indel.tranches -rscriptFile /home/pinedasans/ExomeSeq/VCF/output.indel.plots.R

echo EndRecalibrateIndels
#Apply recalibration from Indels
java -jar GenomeAnalysisTK.jar -T ApplyRecalibration -input /home/pinedasans/ExomeSeq/VCF/joint.output.snp.recal.vcf -o /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.vcf  -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa --ts_filter_level 99.9 -tranchesFile /home/pinedasans/ExomeSeq/VCF/output.indel.tranches -recalFile /home/pinedasans/ExomeSeq/VCF/output.indel.recal -mode INDEL

echo EndApplyRecalibrateIndels
#####


###Quality Control measures

##Exclude the ones that has not passed the filtering
java -jar GenomeAnalysisTK.jar -T SelectVariants -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -V /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.vcf -o /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.vcf --excludeFiltered

##Take only biallelic SNPS
java -jar GenomeAnalysisTK.jar -T SelectVariants -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -V /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.vcf -o /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.biallelic.vcf -restrictAllelesTo BIALLELIC

##Exclude Indels
#java -jar GenomeAnalysisTK.jar -T SelectVariants -R /home/pinedasans/ExomeSeq/BAM/GRCh37.fa -V /home/pinedasans/ExomeSeq/VCF/gVCF/joint.output.snp.indel.recal.filtered.biallelic.vcf -o /home/pinedasans/ExomeSeq/VCF/gVCF/joint.output.snp.indel.recal.filtered.biallelic.SNPs.vcf --selectTypeToExclude INDEL


###To calculate Ts/Tv ratio
export PATH=${PATH}:/home/pinedasans/programs/vcftools_0.1.13/bin/
vcftools --vcf /home/pinedasans/ExomeSeq/VCF/joint.output.snp.indel.recal.filtered.biallelic.vcf --TsTv-summary


