# nf-bowtie2

Create the test directory:
```
mkdir -p ~/nf_atacseq_test/raw_data
```

Download the demo data:
```
cd ~/nf_atacseq_test/raw_data
curl -J -O https://datashare.mpcdf.mpg.de/s/ACJ6T5TVTcvR6fm/download
curl -J -O https://datashare.mpcdf.mpg.de/s/ktjFjaIcLP3lEw0/download

```

Download the paramaters file:
```
cd ~/nf_atacseq_test
curl -J -O https://raw.githubusercontent.com/mpg-age-bioinformatics/nf-bowtie2/main/params.local.json
```

git clone git clone https://github.com/mpg-age-bioinformatics/nf-kallisto.git
git clone https://github.com/mpg-age-bioinformatics/nf-bowtie2.git


Run the workflow:

fastqc
```
PROFILE=raven
nextflow run nf-fastqc -params-file params.local.json -entry images 
nextflow run nf-fastqc -params-file params.local.json
```

flexbar trimming
```
nextflow run nf-flexbar -params-file params.local.json -entry images 
nextflow run nf-flexbar -params-file params.local.json
```

bowtie2
```
nextflow run nf-bowtie2 -params-file params.local.json -entry images  
nextflow run nf-kallisto -params-file params.local.json -entry images  

nextflow run nf-kallisto -params-file  params.local.json -entry get_genome --user "$(id -u):$(id -g)" 

nextflow run nf-kallisto -params-file  params.local.json -entry get_genome  && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry index && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry align && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry mito && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry picard && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry flagstat && \
nextflow run nf-bowtie2 -params-file  params.local.json -entry qccount

```

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```
