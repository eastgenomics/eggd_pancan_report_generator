# pancancer_sample_report

## What does this app do?
Generate aGenerate a html report containing fusion and QC information.

## What are typical use cases for this app?
This app may be executed as a standalone app.

##  What data are required for this app to run?
There are multiple files necessary for this app to run, these can be grouped into dynamic per sample files and static files used for all samples.

Dynamic files:
+ bam
+ fusioninspector abridged
+ rnaseqc metrics
+ rnaseqc coverage
+ rnaseqc exon

Static files:
+ chimeraviz docker
+ CTAT bundle
+ cosmic fusion
+ chimerkb
+ ensdb sqlite
+ capture bed

## What does this app output?
This app outputs a html report.

### This app was made by EMEE GLH