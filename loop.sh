#!/bin/bash

file=/ifs/data/research/projects/kevin/UPD/FamilyInformation_sorted.txt
while read -r line ; do   
  echo line
  echo ${line}
  echo child
  child=$(echo "$line" | cut -f2) 
  echo ${child}
  echo dad
  dad=$(echo "$line" | cut -f3)
  echo ${dad}
  echo mom
  mom=$(echo "$line" | cut -f4)
  echo ${mom}
#  echo complete child
#  echo /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$( echo $line | cut -f2 ).sorted.homREFed.vcf.gz
#  echo complete mom
#  echo /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$( echo $line | cut -f4 ).sorted.homREFed.vcf.gz
  perl -I "/ifs/home/kevin/perl5/lib/perl5/" /ifs/home/kevin/UPDio/version_1.0/UPDio.pl --child_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/${child}.sorted.homREFed.vcf.gz  --mom_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/${mom}.sorted.homREFed.vcf.gz --dad_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/${dad}.sorted.homREFed.vcf.gz --path_to_R /cm/shared/apps/bioinf/R/3.5.1/bin/R --R_scripts_dir /ifs/home/kevin/UPDio/version_1.0/scripts --output_path /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub --include_MI --common_cnv_file /ifs/home/kevin/UPDio/version_1.0/sample_data/common_dels_1percent.tsv  
done < $file
