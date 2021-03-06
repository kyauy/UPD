
# UPD pipeline analysis

## From mendelian inconcordance

### Material and methods

##### Requirements

UPDio
```
R Dependencies
    quantsmooth
    ggplot2
 Perl Dependencies
    Statistics::R (0.31)
    Path::Class
    Vcf
    Iterator::Simple
    List::MoreUtils
    Math::Round
    Const::Fast
UPDio was tested using R version 2.14.1
```

#### With UPDio

One script is missing in version 1.0
```
cp ~/UPDio/version_0.9/scripts/plot_zygosity_and_events.R ~/UPDio/version_1.0/scripts/

```
Pre-processing data

```
for i in /ifs/data/research/projects/kevin/UPD/VCFtest/*.vcf; do cat $i | perl -I "/ifs/home/kevin/perl5/lib/perl5/" ../scripts/sort-vcf | bgzip > $(basename "$i" | cut -d_ -f1).sorted.vcf.gz ; perl -I "/ifs/home/kevin/perl5/lib/perl5/" ../scripts/add_hom_refs_to_vcf.pl --polymorphic_sites ../sample_data/common_variants_within_well_covered_target_regions.txt --no_homREF_vcf $(basename "$i" | cut -d_ -f1).sorted.vcf.gz | bgzip > $(basename "$i" | cut -d_ -f1).sorted.homREFed.vcf.gz ; done
```

Get CNV
```
for i in /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/CNVCalling*_conifer/*_calls.txt ; do cp -n $i /ifs/data/research/projects/kevin/UPD/CNV ; done
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA*/GRCh37_*/CNVCalling*_conifer/*_calls.txt ; do cp -n $i /ifs/data/research/projects/kevin/UPD/CNV ; done

### Select sample to use CNV

for i in $(cat /ifs/data/research/projects/kevin/UPD/results/compare/upd_called.txt) ; do cp -n /ifs/data/research/projects/kevin/UPD/CNV/$(echo $i)_calls.txt /ifs/data/research/projects/kevin/UPD/CNV_pub ; done
for i in $(cat /ifs/data/research/projects/kevin/UPD/results/compare/upd_called.txt) ; do grep $i /ifs/data/research/projects/kevin/UPD/FamilyInformation_sorted.txt >> /ifs/data/research/projects/kevin/UPD/FamilyInformation_sorted_UPD_wo_CNV.txt ; done

for i in /ifs/data/research/projects/kevin/UPD/CNV_pub/*calls.txt ; do tail -n+2 $i | cut -f2,3,4 > /ifs/data/research/projects/kevin/UPD/CNV_pub/$(basename $i | cut -d_ -f1).cnv.txt ; done
```

Run UPDio

```
perl -I "/ifs/home/kevin/perl5/lib/perl5/" UPDio.pl --child_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA18-03831B.sorted.homREFed.vcf.gz --mom_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA18-03833B.sorted.homREFed.vcf.gz --dad_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA18-03832B.sorted.homREFed.vcf.gz --path_to_R /cm/shared/apps/bioinf/R/3.5.1/bin/R --R_scripts_dir /ifs/home/kevin/UPDio/version_1.0/scripts --include_MI --common_cnv_file /ifs/home/kevin/UPDio/version_1.0/sample_data/common_dels_1percent.tsv --increase_cnv_filtering --child_cnv_data /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA18-03831B_cnv.txt

perl -I "/ifs/home/kevin/perl5/lib/perl5/" UPDio.pl --child_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA11-26540B.sorted.homREFed.vcf.gz --mom_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA12-07537B.sorted.homREFed.vcf.gz --dad_vcf /ifs/home/kevin/UPDio/version_1.0/pre_processing/DNA12-07536B.sorted.homREFed.vcf.gz --path_to_R /cm/shared/apps/bioinf/R/3.5.1/bin/R --R_scripts_dir /ifs/home/kevin/UPDio/version_1.0/scripts
```

Two loops that rule them all

```
## PREPROCESSING LOOP
for i in /ifs/data/research/projects/kevin/UPD/VCF/*.vcf; do cat $i | perl -I "/ifs/home/kevin/perl5/lib/perl5/" /ifs/home/kevin/UPDio/version_1.0/scripts/sort-vcf | bgzip > /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(basename "$i" | cut -d_ -f1).sorted.vcf.gz ; perl -I "/ifs/home/kevin/perl5/lib/perl5/"  /ifs/home/kevin/UPDio/version_1.0/scripts/add_hom_refs_to_vcf.pl --polymorphic_sites  /ifs/home/kevin/UPDio/version_1.0/sample_data/common_variants_within_well_covered_target_regions.txt --no_homREF_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(basename "$i" | cut -d_ -f1).sorted.vcf.gz | bgzip > /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(basename "$i" | cut -d_ -f1).sorted.homREFed.vcf.gz ; rm -f /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(basename "$i" | cut -d_ -f1).sorted.vcf.gz ; done

## UPDio LOOP
for i in $(cat < "/ifs/data/research/projects/kevin/UPD/family_sorted_mendel.txt") ; do IFS=$'\n' ; perl -I "/ifs/home/kevin/perl5/lib/perl5/" UPDio.pl --child_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f2).sorted.homREFed.vcf.gz  --mom_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f4).sorted.homREFed.vcf.gz --dad_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f3).sorted.homREFed.vcf.gz --path_to_R /cm/shared/apps/bioinf/R/3.5.1/bin/R --R_scripts_dir /ifs/home/kevin/UPDio/version_1.0/scripts --output_path ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed --include_MI --common_cnv_file /ifs/home/kevin/UPDio/version_1.0/sample_data/common_dels_1percent.tsv ; done
```

Plot legend

```
BPI = Biparental
UA = Uniparental Ambigous
UI = Uniparental Isodisomic

Informative genotype combinations. Sites at which parents are opposing homozygotes
and the child is heterozygous are diagnostic of biparental inheritance. Uniparental
inheritance combinations include those that obligately result from isodisomy (UI), and
those that may result from either heterodisomy or isodisomy (UA) as the proband
alleles may have arisen from a duplication of one parental homolog, or may present
both homologs. True isodisomy events will have mixtures of both of these informative
types, while true heterodisomy events will be void of UI types.
```

Create a file with all data

```
for i in /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed/*upd; do awk '{print FILENAME (NF?"\t":"") $0}' $i | sed 's/.sorted.homREFed.upd//g' >> /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed/complete_UPD.tsv ; done
```

Specific analysis with trio for publication
```
sed 's/\t.\{5\}/&-/g' FamilyInformation.txt | sed '/^ /d' | sed '/^?/d' | sort -n > FamilyInformation_sorted.txt
for i in $(cat < "/ifs/data/research/projects/kevin/UPD/FamilyInformation_sorted.txt") ; do IFS=$'\n' ; perl -I "/ifs/home/kevin/perl5/lib/perl5/" /ifs/home/kevin/UPDio/version_1.0/UPDio.pl --child_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f2).sorted.homREFed.vcf.gz  --mom_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f4).sorted.homREFed.vcf.gz --dad_vcf /ifs/data/research/projects/kevin/UPD/VCF_UPDio/$(echo $i | cut -f3).sorted.homREFed.vcf.gz --path_to_R /cm/shared/apps/bioinf/R/3.5.1/bin/R --R_scripts_dir /ifs/home/kevin/UPDio/version_1.0/scripts --output_path ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub --include_MI --common_cnv_file /ifs/home/kevin/UPDio/version_1.0/sample_data/common_dels_1percent.tsv ; done
## UPDio loop
sh loop.sh
## UPDio loop for UPD with CNV data
sh loop_w_cnv.sh

## file with all data
for i in /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub/*upd; do awk '{print FILENAME (NF?"\t":"") $0}' $i | sed 's/.sorted.homREFed.upd//g' | sed 's/\/ifs\/data\/research\/projects\/kevin\/UPD\/VCF_UPDio_processed_pub\///g' >> /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub/complete_UPD.tsv ; done

## file with all data + cnv
for i in /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub_cnv/*upd; do awk '{print FILENAME (NF?"\t":"") $0}' $i | sed 's/.sorted.homREFed.upd//g' | sed 's/\/ifs\/data\/research\/projects\/kevin\/UPD\/VCF_UPDio_processed_pub_cnv\///g' >> /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub_cnv/complete_UPD_pub_cnv.tsv ; done

## COPY png
for i in /ifs/data/research/projects/kevin/UPD/VCF_UPDio_processed_pub/*png; do cp -n $i /ifs/data/research/projects/kevin/UPD/results_png ; done
## Compare with ROH
for i in $(cat results/ROH_complete_iso.txt) ; do ls VCF_UPDio_processed_pub | grep $i ; done
```

Selected UPD from 4912 trio
```
### list of all samples
15:43:28 kevin::login03 { ~/UPD/data/Whole/VCF_UPDio_processed_cnv }-> ls *.upd | cut -f1 -d. > allsamples.txt

### list of selected upd
11:34:25 kevin::login03 { ~/UPD/data/Whole/results_png }-> ls | grep -f selectedupd.txt > selectedupd_filename.txt
11:34:51 kevin::login03 { ~/UPD/data/Whole/results_png }-> mkdir selected
11:35:10 kevin::login03 { ~/UPD/data/Whole/results_png }-> cat selectedupd_filename.txt |  xargs -I % cp % selected/

### get H3M2 plot of selected upd
while read -r line; do
    cp -n /ifs/data/diagnostics/bgi/exomes/work/${line}/GRCh37_*/ROH_*/*.pdf selected/
    cp -n /ifs/data/diagnostics/nextseq/exomes/work/${line}/GRCh37_*/ROH_*/*.pdf selected/
done < selectedupd.txt


### get MAD score only for UPDIO trio
15:49:46 kevin::login03 { ~/UPD/data/Whole/VCF_UPDio_processed_cnv }-> head -n1 ../results/allROH_processed.tsv > trioupdioROH_processed_MADscore.tsv
15:51:07 kevin::login03 { ~/UPD/data/Whole/VCF_UPDio_processed_cnv }-> grep -f ~/UPD/data/Whole/VCF_UPDio_processed_cnv/allsamples.txt ../results/allROH_processed_MADscore_quality.tsv >> trioupdioROH_processed_MADscore.tsv
```

#### Complete process with bcftools

Prepare data for mendelian concordance analysis
```
### Prepare file with child and parents name
16:06:46 kevin::login03 { ~/UPD/data/Whole }-> sed "s/.\{5\}/&-/" family.txt | sed 's/\t.\{5\}/&-/g' | sed '/^ /d' | sed '/^?/d' | sed '1d' | sort | sed '1d' > family_sorted.txt

### Fit for PED format in excel and remove duplicate: family_sorted_mendel.txt

### Remove double -- : sed -i 's/--/-/g' family_sorted_mendel.txt

### Remove duplicate
17:38:36 kevin::login01 { ~/UPD/data/Whole }-> uniq /ifs/data/research/projects/kevin/UPD/child_sorted.txt >  /ifs/data/research/projects/kevin/UPD/child_sorted_uniq.txt

### Create PED for each family
for i in $(cat < "/ifs/data/research/projects/kevin/UPD/child_sorted_uniq.txt") ; do grep "$i" /ifs/data/research/projects/kevin/UPD/family_sorted_mendel.txt > /ifs/data/research/projects/kevin/UPD/VCFtools/$(echo "$i").mendel; done

### Copy all needed VCF

for i in /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/SNVC*/*.vcf* ; do cp -n $i /ifs/data/research/projects/kevin/UPD/VCF ; done
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA*/GRCh37_*/SNVC*/*.vcf* ; do cp -n $i /ifs/data/research/projects/kevin/UPD/VCF ; done


```

Merge VCF

```
### Create new folder with fitted for merge VCF

cp -r /ifs/data/research/projects/kevin/UPD/VCF /ifs/data/research/projects/kevin/UPD/VCF_MC

for file in *haplotypecaller*; do mv "$file" "${file/_haplotypecaller/}" ; done
for file in *unifiedgenotyper*; do mv "$file" "${file/_unifiedgenotyper/}" ; done

## Add - for DNA names in VCF
for i in /ifs/data/research/projects/kevin/UPD/VCF_MC ; do sed -i "s/FORMAT\t.\{5\}/&-/" $i ; sed -i 's/--/-/g' $i ; done

### Merge in VCFmerge folder
mkdir /ifs/data/research/projects/kevin/UPD/VCFmerge/
for i in $(cat < "/ifs/data/research/projects/kevin/UPD/family_sorted_mendel.txt") ; do IFS=$'\n' ; bcftools merge -o /ifs/data/research/projects/kevin/UPD/VCFmerge/$(echo $i | cut -f2)_merged.vcf -O v -m none /ifs/data/research/projects/kevin/UPD/VCF_MC/$(echo $i | cut -f2).vcf /ifs/data/research/projects/kevin/UPD/VCF_MC/$(echo $i | cut -f3).vcf /ifs/data/research/projects/kevin/UPD/VCF_MC/$(echo $i | cut -f4).vcf ; done
```

#### Collect data from previous pipeline

Copy all data

```
for i in /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/SNVC*/SNVA*/SNVD*/*mendelian_inheritance.txt ; do cp -n $i /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw ; done

for i in /ifs/data/diagnostics/bgi/exomes/work/DNA*/GRCh37_*/SNVC*/SNVA*/SNVD*/*mendelian_inheritance.txt ; do cp -n $i /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw ; done
```

Get processable data for z-score.
Beware : I had a one issue from DNA17-06937B_mendelian_inheritance.txt as no variant in chrX were inherited from the father.

```
### Remove chrY and chrM since they have 0 and we can't divide by zero
for i in /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw/*_mendelian_inheritance.txt; do sed -i '/chrY/d' $i ; sed -i '/chrM/d' $i ; done

### Remove header line
for i in /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw/*_mendelian_inheritance.txt; do sed -i '1d' $i ; done

### Get data.frame for maternal UPD
for i in /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw/*mendelian_inheritance.txt  ; do printf "Sample\n$(basename "$i" | cut -d_ -f1)" > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3/$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

### Get data.frame for paternal UPD
for i in /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_raw/*mendelian_inheritance.txt  ; do printf "Sample\n$(basename "$i" | cut -d_ -f1)" > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($2/$3)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done
```


Create complete dataset with all patients
```
### maternal

head -n1 /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/DNA18-16505B_mendelian_inheritance.processed.bed  > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/all_mendelian_inheritance_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> all_mendelian_inheritance_processed.tsv

### paternal

head -n1 /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/DNA18-16505B_mendelian_inheritance.processed.bed  > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/all_mendelian_inheritance_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> all_mendelian_inheritance_processed.tsv


```

#### Create z-score for all chromosome

Create z-score and column for filtering
```
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/all_mendelian_inheritance_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/all_mendelian_inheritance_processed.tsv

```

Add column for filtering

```
awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/all_mendelian_inheritance_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_processed/all_mendelian_inheritance_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/all_mendelian_inheritance_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/mendelian_inheritance_paternal_processed/all_mendelian_inheritance_paternal_processed_zscore_quality.tsv

```

Filter file with =< X chromosome under 2
```
### ALL
cut -f1-24 all_mendelian_inheritance_paternal_processed_zscore_quality.tsv > all_mendelian_inheritance_paternal_processed_zscore_quality_R.tsv
cut -f1-24 all_mendelian_inheritance_maternal_processed_zscore_quality.tsv > all_mendelian_inheritance_maternal_processed_zscore_quality_R.tsv


### P1
head -n1 all_mendelian_inheritance_paternal_processed_zscore_quality.tsv > all_mendelian_inheritance_paternal_processed_zscore_quality_1.tsv
tail -n+2 all_mendelian_inheritance_paternal_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) == 1 {print $0}' >> all_mendelian_inheritance_paternal_processed_zscore_quality_1.tsv
cut -f1-24 all_mendelian_inheritance_paternal_processed_zscore_quality_1.tsv > all_mendelian_inheritance_paternal_processed_zscore_quality_1_R.tsv

#### M1
head -n1 all_mendelian_inheritance_maternal_processed_zscore_quality.tsv > all_mendelian_inheritance_maternal_processed_zscore_quality_1.tsv
tail -n+2 all_mendelian_inheritance_maternal_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) == 1 {print $0}' >> all_mendelian_inheritance_maternal_processed_zscore_quality_1.tsv
cut -f1-24 all_mendelian_inheritance_maternal_processed_zscore_quality_1.tsv > all_mendelian_inheritance_maternal_processed_zscore_quality_1_R.tsv

### P2
head -n1 all_mendelian_inheritance_paternal_processed_zscore_quality.tsv > all_mendelian_inheritance_paternal_processed_zscore_quality_2.tsv
tail -n+2 all_mendelian_inheritance_paternal_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) <= 2 && $(NF-3) > 0 {print $0}' >> all_mendelian_inheritance_paternal_processed_zscore_quality_2.tsv
cut -f1-24 all_mendelian_inheritance_paternal_processed_zscore_quality_2.tsv > all_mendelian_inheritance_paternal_processed_zscore_quality_2_R.tsv

#### M2
head -n1 all_mendelian_inheritance_maternal_processed_zscore_quality.tsv > all_mendelian_inheritance_maternal_processed_zscore_quality_2.tsv
tail -n+2 all_mendelian_inheritance_maternal_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) <= 2 && $(NF-3) > 0 {print $0}' >> all_mendelian_inheritance_maternal_processed_zscore_quality_2.tsv
cut -f1-24 all_mendelian_inheritance_maternal_processed_zscore_quality_2.tsv > all_mendelian_inheritance_maternal_processed_zscore_quality_2_R.tsv

```


## From ROH H3M2 WES data

### Material and methods

#### Requirements

R

#### Collect data

Number of samples :
```
15:59:10 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA13-*/GRCh37_*/ROH_*/*.bed | wc -l
2741
16:00:10 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA12-*/GRCh37_*/ROH_*/*.bed | wc -l
904
16:00:49 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA11-*/GRCh37_*/ROH_*/*.bed | wc -l
382
16:01:02 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA10-*/GRCh37_*/ROH_*/*.bed | wc -l
253
16:01:12 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA09-*/GRCh37_*/ROH_*/*.bed | wc -l
178
16:01:21 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA08-*/GRCh37_*/ROH_*/*.bed | wc -l
115
16:01:28 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA07-*/GRCh37_*/ROH_*/*.bed | wc -l
55
16:01:33 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA06-*/GRCh37_*/ROH_*/*.bed | wc -l
63
16:01:38 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA05-*/GRCh37_*/ROH_*/*.bed | wc -l
45
16:01:43 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA03-*/GRCh37_*/ROH_*/*.bed | wc -l
16
16:01:47 kevin::login03 { ~ }-> ls -f /ifs/data/diagnostics/bgi/exomes/work/DNA02-*/GRCh37_*/ROH_*/*.bed | wc -l
7
15:40:35 kevin::login02 { /ifs/data/research/projects/kevin/UPD/results }-> ls  /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/ROH_*/*.bed | wc -l
594
```

Get name for proband
```
for i in /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/*.properties; do grep "family_fraction_numbers" $i >>  /ifs/data/research/projects/kevin/UPD/samples.txt ; done
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA*/GRCh37_*/*.properties; do grep "family_fraction_numbers" $i >>  /ifs/data/research/projects/kevin/UPD/samples.txt ; done

17:05:51 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sort samples.txt | uniq | sed 's/family_fraction_numbers=//g' | sed 's/;/\t/g' | sed 's/\\//g' | sed 's/KIND://g' | sed 's/VADER://g' | sed 's/MOEDER://g' >> family.txt

#### split
tail -n +2 /ifs/data/research/projects/kevin/UPD/family.txt | cut -f1 > /ifs/data/research/projects/kevin/UPD/child.txt
tail -n +2 /ifs/data/research/projects/kevin/UPD/family.txt | cut -f2 > /ifs/data/research/projects/kevin/UPD/father.txt
tail -n +2 /ifs/data/research/projects/kevin/UPD/family.txt | cut -f3 > /ifs/data/research/projects/kevin/UPD/mother.txt

17:31:31 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sed -i -e "s/.\{5\}/&-/" child.txt
17:32:05 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sed -i -e "s/.\{5\}/&-/" father.txt
17:32:18 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sed -i -e "s/.\{5\}/&-/" mother.txt

17:36:39 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sort child.txt | sed '/^ /d' | sed '/^?/d' | sed '1d' > child_sorted.txt
17:39:33 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sort father.txt | sed '/^ /d' | sed '/^?/d' | sed '1d' > father_sorted.txt
17:39:41 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> sort mother.txt | sed '/^ /d' | sed '/^?/d' | sed '1d' > mother_sorted.txt

```


Copy all ROH file done from known data
```

###### Home

for i in /ifs/data/diagnostics/nextseq/exomes/work/DNA*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/nextseq_raw/ ; done

###### BGI
## before 2012
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA0*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/before2012_ROH_raw/ ; done
DNA0*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA10-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/before2012_ROH_raw/ ; done
DNA10*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA11-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/before2012_ROH_raw/ ; done
DNA1*/GRCh37_*/ROH_*/*.bed

## after 2012
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA12-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2012_ROH_raw/ ; done
DNA12*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA13-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2013_ROH_raw/ ; done
DNA13*/GRCh37_*/ROH_*/*.bed

for i in /ifs/data/diagnostics/bgi/exomes/work/DNA14-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2014_ROH_raw/ ; done
DNA14*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA15-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2015_ROH_raw/ ; done
DNA15*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA16-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2016_ROH_raw/ ; done
DNA16*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA17-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2017_ROH_raw/ ; done
DNA17*/GRCh37_*/ROH_*/*.bed
for i in /ifs/data/diagnostics/bgi/exomes/work/DNA18-*/GRCh37_*/ROH_*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/2018_ROH_raw_all/ ; done
DNA18*/GRCh37_*/ROH_*/*.bed

#### ALL
for i in /ifs/data/research/projects/kevin/UPD/*ROH_raw*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/ROH_raw_all/ ; done
for i in /ifs/data/research/projects/kevin/UPD/*nextseq_raw*/*.bed ; do cp -n $i /ifs/data/research/projects/kevin/UPD/ROH_raw_all/ ; done


```

Get length by chromosome and processable data for z-score calling
```

### Nextseq

for i in /ifs/data/research/projects/kevin/UPD/nextseq_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

#### BGI
## Before 2012

for i in /ifs/data/research/projects/kevin/UPD/before2012_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

### Since 2012

for i in /ifs/data/research/projects/kevin/UPD/2012_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2013_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done


for i in /ifs/data/research/projects/kevin/UPD/2014_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2015_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2016_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2017_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2018_ROH_raw_all/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

### ALL

for i in /ifs/data/research/projects/kevin/UPD/ROH_raw_all/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done
```

Create whole dataset for all years
```
#### NextSeq

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/nextseq_processed/nextseq_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/nextseq_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> nextseq_processed.tsv

#### BGI
## Before 2012
head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/mendelian_inheritance_processed/all_mendelian_inheritance_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> before2012_allROH_processed.tsv

## After 2012

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/2012_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2012_allROH_processed.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/2013_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2013_allROH_processed.tsv


head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/2014_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2014_allROH_processed.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/2015_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2015_allROH_processed.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/2016_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2016_allROH_processed.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/2017_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2017_allROH_processed.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/2018_allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> 2018_allROH_processed.tsv

### ALL

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/ROH_all_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> allROH_processed.tsv

```

Jakob tricks and advices :
Use SNP, get the biggest and compare.
```
#### use SNP and get only big ones
11:20:06 kevin::login02 { /ifs/data/research/projects/kevin/UPD/ROH_raw_all }-> for i in *.bed; do cut -f5 $i >> allsnp_number.txt ; done
11:24:29 kevin::login02 { /ifs/data/research/projects/kevin/UPD/ROH_raw_all }-> cat allsnp_number.txt | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
       0.0   195.0   259.0   330.7   361.0 36420.0
11:24:50 kevin::login02 { /ifs/data/research/projects/kevin/UPD/ROH_raw_all }-> cat allsnp_number.txt | Rscript -e 'quantile (as.numeric (readLines ("stdin")), probs=c(0.025, 0.5, 0.975))'
 2.5%   50% 97.5%
   123   259   859

### get for all
for i in /ifs/data/research/projects/kevin/UPD/ROH_raw_all/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '$5 > 361 {l[$1]+=($5)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/$(basename "$i" | cut -d. -f1).sample.bed ; done

### Merge
head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/allROH_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> allROH_processed.tsv

## R
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/allROH_processed.tsv

## filter

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/ROH_all_processed_jg/allROH_processed_jg_zscore_quality.tsv

### ALL
cut -f1-24 allROH_processed_jg_zscore_quality.tsv > allROH_processed_jg_zscore_quality_R.tsv

### Filter
#### 2
head -n1 allROH_processed_zscore_quality.tsv > allROH_processed_zscore_quality_filtered_2.tsv
tail -n+2 allROH_processed_jg_zscore_quality.tsv | awk -F "\t" '$(NF-3) <= 2 && $(NF-3) > 0 {print $0}' >> allROH_processed_jg_zscore_quality_filtered_2.tsv
cut -f1-24 allROH_processed_jg_zscore_quality_filtered_2.tsv > allROH_processed_jg_zscore_quality_filtered_2_R.tsv
```

Get addition size of all ROH by sample

```
15:56:01 kevin::login03 { ~/UPD/data/Whole/ROH_all_processed }-> awk '{for(i=1;i<=NF;i++) t+=$i; print $0"\t"t; t=0}' allROH_processed.tsv | cut -f1,26 | sed 's/e\t0/e\tSize/' > allROH_processed_sumsize.tsv
```

#### Create z-score for all chromosome

Create z-score and column for filtering
```
### NextSeq
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/nextseq_processed/nextseq_processed.tsv

### BGI
## Before 2012
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/before2012_allROH_processed.tsv

## After 2012
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/2012_allROH_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/2013_allROH_processed.tsv

Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/2014_allROH_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/2015_allROH_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/2016_allROH_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/2017_allROH_processed.tsv
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/2018_allROH_processed.tsv

### ALL
Rscript ~/UPD/zscore.R /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed.tsv

```
Add column for filtering

```
### NextSeq

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/nextseq_processed/nextseq_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/nextseq_processed/nextseq_processed_zscore_quality.tsv

### BGI
## Before 2012

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/before2012_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/before2012_allROH_processed_zscore_quality.tsv


#### After 2012

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/2012_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/2012_allROH_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/2013_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/2013_allROH_processed_zscore_quality.tsv



awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/2014_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/2014_allROH_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/2015_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/2015_allROH_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/2016_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/2016_allROH_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/2017_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/2017_allROH_processed_zscore_quality.tsv

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/2018_allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/2018_allROH_processed_zscore_quality.tsv

### ALL

awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv

15:40:00 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> for i in */*quality.tsv ; do cp $i results/ ; done

```
Split for child, father and mother sample

```
tail -n +2 /ifs/data/research/projects/kevin/UPD/results/allROH_processed_zscore_quality.tsv | sort > /ifs/data/research/projects/kevin/UPD/results/allROH_processed_zscore_quality_sorted.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv > /ifs/data/research/projects/kevin/UPD/results/child_allROH_processed_zscore_quality.tsv
grep -f /ifs/data/research/projects/kevin/UPD/child_sorted.txt /ifs/data/research/projects/kevin/UPD/results/allROH_processed_zscore_quality_sorted.tsv  >> /ifs/data/research/projects/kevin/UPD/results/child_allROH_processed_zscore_quality.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv > /ifs/data/research/projects/kevin/UPD/results/father_allROH_processed_zscore_quality.tsv
grep -f /ifs/data/research/projects/kevin/UPD/father_sorted.txt /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv >> /ifs/data/research/projects/kevin/UPD/results/father_allROH_processed_zscore_quality.tsv

head -n1 /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv > /ifs/data/research/projects/kevin/UPD/results/mother_allROH_processed_zscore_quality.tsv
grep -f /ifs/data/research/projects/kevin/UPD/mother_sorted.txt /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_zscore_quality.tsv >> /ifs/data/research/projects/kevin/UPD/results/mother_allROH_processed_zscore_quality.tsv
```

Filter file with =< X chromosome under 2
```
### ALL
cut -f1-24 allROH_processed_zscore_quality.tsv > allROH_processed_zscore_quality_R.tsv

### >2 - 1
head -n1 allROH_processed_zscore_quality.tsv > allROH_processed_zscore_quality_filtered_1.tsv
tail -n+2 allROH_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) == 1 {print $0}' >> allROH_processed_zscore_quality_filtered_1.tsv
cut -f1-24 allROH_processed_zscore_quality_filtered_1.tsv > allROH_processed_zscore_quality_filtered_1_R.tsv

#### >2 - 2
head -n1 allROH_processed_zscore_quality.tsv > allROH_processed_zscore_quality_filtered_2.tsv
tail -n+2 allROH_processed_zscore_quality.tsv | awk -F "\t" '$(NF-3) <= 2 && $(NF-3) > 0 {print $0}' >> allROH_processed_zscore_quality_filtered_2.tsv
cut -f1-24 allROH_processed_zscore_quality_filtered_2.tsv > allROH_processed_zscore_quality_filtered_2_R.tsv

### > 1.5 and 2- 1
head -n1 allROH_processed_zscore_quality.tsv > allROH_processed_zscore_quality_filtered_1.5.tsv
tail -n+2 allROH_processed_zscore_quality.tsv | awk -F "\t" '$(NF-4) == 1 {print $0}' | awk -F "\t" '$(NF-3) == 1 {print $0}'  >> allROH_processed_zscore_quality_filtered_1.5.tsv
cut -f1-24 allROH_processed_zscore_quality_filtered_1.5.tsv > allROH_processed_zscore_quality_filtered_1.5_R.tsv


```

#### Create median absolute deviation for all chromosome

Create MAD and column for filtering
Need the quantable package
```
### ALL
Rscript ~/UPD/MADscore.R /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed.tsv

### Create columns for filtering
# without median awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-1;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-2;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-3;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-4;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore_quality.tsv

## ALL samples
awk -F "\t" '{if(NR==1) c="Number>1.5"; else for(i=2;i<=NF-1;i++) c+=($i>1.5); print $0,c; c=0}' /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore.tsv | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>2"; else for(i=2;i<=NF-2;i++) d+=($i>2); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>3"; else for(i=2;i<=NF-3;i++) d+=($i>3); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>4"; else for(i=2;i<=NF-4;i++) d+=($i>4); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | awk -F "\t" '{if(NR==1) d="Number>5"; else for(i=2;i<=NF-5;i++) d+=($i>5); print $0,d; d=0}' | column -t | sed 's/ \+ /\t/g' | sed '1s/^/Sample\t/'  > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore_quality.tsv
cp /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore_quality.tsv /ifs/data/research/projects/kevin/UPD/results

## Just publication trio
head -n1 /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore_quality.tsv > /ifs/data/research/projects/kevin/UPD/ROH_all_processed/pubtrioROH_processed_MADscore_quality.tsv
for i in $( cut -f2 /ifs/data/research/projects/kevin/UPD/FamilyInformation_sorted.txt ); do grep $i /ifs/data/research/projects/kevin/UPD/ROH_all_processed/allROH_processed_MADscore_quality.tsv >> /ifs/data/research/projects/kevin/UPD/ROH_all_processed/pubtrioROH_processed_MADscore_quality.tsv ; done

### Filter file with =< X chromosome under 2
### ALL
cut -f1-24 allROH_processed_MADscore_quality.tsv > allROH_processed_MADscore_quality_R.tsv

#### >2 - 2
head -n1 allROH_processed_MADscore_quality.tsv > allROH_processed_MADscore_quality_filtered_2.tsv
tail -n+2 allROH_processed_MADscore_quality.tsv | awk -F "\t" '$(NF-3) <= 2 && $(NF-3) > 0 {print $0}' >> allROH_processed_MADscore_quality_filtered_2.tsv
cut -f1-24 allROH_processed_MADscore_quality_filtered_2.tsv > allROH_processed_MADscore_quality_filtered_2_R.tsv

#### >3 - 2
head -n1 allROH_processed_MADscore_quality.tsv > allROH_processed_MADscore_quality_filtered_3.tsv
tail -n+2 allROH_processed_MADscore_quality.tsv | awk -F "\t" '$(NF-2) <= 2 && $(NF-2) > 0 {print $0}' >> allROH_processed_MADscore_quality_filtered_3.tsv
cut -f1-23 allROH_processed_MADscore_quality_filtered_3.tsv > allROH_processed_MADscore_quality_filtered_3_R.tsv

### for publication
head -n1 pubtrioROH_processed_MADscore_quality.tsv > pubtrioROH_processed_MADscore_quality_filtered_3.tsv
tail -n+2 pubtrioROH_processed_MADscore_quality.tsv | awk -F "\t" '$(NF-2) <= 2 && $(NF-2) > 0 {print $0}' >> pubtrioROH_processed_MADscore_quality_filtered_3.tsv
cut -f1-23 pubtrioROH_processed_MADscore_quality_filtered_3.tsv > pubtrioROH_processed_MADscore_quality_filtered_3_R.tsv

```

## Compare ROH vs UPDio

Good ROH on UPDio
```
for i in $( cut -f2 /ifs/data/research/projects/kevin/UPD/results/compare/pubtrioROH_med3_2.list ); do grep $i /ifs/data/research/projects/kevin/UPD/results/complete_UPD_pub.tsv >> /ifs/data/research/projects/kevin/UPD/results/compare/pubtrioROH_med3_2.updio.tsv ; done
```

Good UPDio Isodisomy on ROH
```
sed -i 's/ \+ /\t/g' complete_UPD_pub_cnv_formerge.tsv
sed -i 's/ /\t/g' complete_UPD_pub_cnv_formerge.tsv
19:14:47 kevin::login04 { ~/UPD/data/Whole/results/compare }-> python merge.py
```

