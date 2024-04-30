# nf-bowtie2

Create the test directory:
```
mkdir -p /tmp/nextflow_atac_local_test/raw_data
```

Download the demo data:
```
cd /tmp/nextflow_atac_local_test/raw_data
curl -J -O https://datashare.mpcdf.mpg.de/s/ACJ6T5TVTcvR6fm/download
curl -J -O https://datashare.mpcdf.mpg.de/s/ktjFjaIcLP3lEw0/download

```

Download the paramaters file:
```
cd /tmp/nextflow_atac_local_test
PARAMS=params.local.json
curl -J -O https://raw.githubusercontent.com/mpg-age-bioinformatics/nf-bowtie2/main/${PARAMS}
```

Get the latest repo:
```
git clone https://github.com/mpg-age-bioinformatics/nf-fastqc.git
git clone https://github.com/mpg-age-bioinformatics/nf-kallisto.git
git clone https://github.com/mpg-age-bioinformatics/nf-bowtie2.git
```

Run the workflow:

fastqc
```
nextflow run nf-fastqc -params-file ${PARAMS} -entry images --user "$(id -u):$(id -g)"  
nextflow run nf-fastqc -params-file ${PARAMS} --user "$(id -u):$(id -g)"  
```

flexbar trimming
```
nextflow run nf-flexbar -params-file ${PARAMS} -entry images --user "$(id -u):$(id -g)"  
nextflow run nf-flexbar -params-file ${PARAMS} --user "$(id -u):$(id -g)"  
```

bowtie2
```
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry images --user "$(id -u):$(id -g)"
nextflow run nf-kallisto -params-file ${PARAMS} -entry images --user "$(id -u):$(id -g)"

nextflow run nf-kallisto -params-file ${PARAMS} -entry get_genome --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry index --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry align --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry mito --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry picard --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry flagstat --user "$(id -u):$(id -g)" && \
nextflow run nf-bowtie2 -params-file ${PARAMS} -entry qccount --user "$(id -u):$(id -g)" 
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
