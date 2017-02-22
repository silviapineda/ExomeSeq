with open('/home/pinedasans/ExomeSeq/1000G/variantList_vcf.txt') as f:
    lines1 = f.readlines()
    #print (lines1)
  
fout = open('/home/pinedasans/ExomeSeq/1000G/1000G_genotypes_variantList.vcf','w')

for i in range(1,22):
    #print (i)
    with open("/home/pinedasans/ExomeSeq/1000G/ALL.chr" + str(i) + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf") as f:
        print (f)
        for line in f:
            if any(line.startswith(x.strip()) for x in lines1):
                #print (line)
                fout.write( line + '\n')
        
        