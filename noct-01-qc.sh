#!/usr/env/bin bash

source noct-00-envSetup.sh

mkdir "$qcDir"

mkdir "$qcDir"fastqc/
fastqc "$seqDir"*"$fastqSuffix" -o "$qcDir"fastqc/ --noextract
mkdir "$qcDir"multiqc/
multiqc "$qcDir"fastqc/* -o "$qcDir"multiqc/
