with open('/~/ExomeSeq/1000G/variantList_vcf.txt') as f:
    lines1 = f.readlines()
    #print lines1
  
fout = open('/~/ExomeSeq/1000G/1000G_genotypes_variantList.vcf','w')

with open('/~/ExomeSeq/1000G/1000G_genotypes.vcf') as f:
  for line in f:
    if any(line.startswith(x.strip()) for x in lines1):
        print line
        fout.write( line + '\n')