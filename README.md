
# UPD

## from mendelian inconcordance

### Material and methods

#### Requirements

#### Collect data for reference panel



## from ROH H3M2 WES data

### Material and methods

#### Requirements

R

#### Collect data for reference panel

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

```

Get length by chromosome and processable data for z-score calling
```

### Nextseq

for i in /ifs/data/research/projects/kevin/UPD/nextseq_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/nextseq_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

#### BGI
### Before 2012

for i in /ifs/data/research/projects/kevin/UPD/before2012_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

### Since 2012

for i in /ifs/data/research/projects/kevin/UPD/2012_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2012_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2013_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2013_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done


for i in /ifs/data/research/projects/kevin/UPD/2014_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2014_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2015_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2015_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2016_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2016_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2017_ROH_raw/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2017_ROH_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done

for i in /ifs/data/research/projects/kevin/UPD/2018_ROH_raw_all/*.bed ; do printf "Sample\n$(basename "$i" | cut -d. -f1)" > /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; awk -F "\t" '{l[$1]+=($3-$2)} END {for (i in l) {print i,l[i]}}' "$i" | sort -k1,1V | sed 's/ /\t/g' | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | paste /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed - > /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).processed.bed ; rm -f /ifs/data/research/projects/kevin/UPD/2018_ROH_all_processed/$(basename "$i" | cut -d. -f1).sample.bed ; done
```

Create whole dataset for all years
```
#### NextSeq

head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/nextseq_processed/nextseq_processed.tsv
cd /ifs/data/research/projects/kevin/UPD/nextseq_processed/
find . -name "*.bed" | xargs -n 1 tail -n +2 >> nextseq_processed.tsv

#### BGI
## Before 2012
head -n1 /ifs/data/research/projects/kevin/UPD/2018_ROH_processed/DNA18-16523B.processed.bed >  /ifs/data/research/projects/kevin/UPD/before2012_ROH_processed/before2012_allROH_processed.tsv
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

15:40:00 kevin::login02 { /ifs/data/research/projects/kevin/UPD }-> for i in */*quality.tsv ; do cp $i results/ ; done

```
