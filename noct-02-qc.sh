#!/usr/env/bin bash

source noct-00-envSetup.sh

# mkdir "$qcDir"

qcStage="cutAndTrim_"

mkdir "$qcDir""$qcStage"fastqc/
fastqc "$cutDir"*.paired"$trimSuffix" -o "$qcDir""$qcStage"fastqc/ --noextract
mkdir "$qcDir""$qcStage"multiqc/
multiqc "$qcDir""$qcStage"fastqc/* -o "$qcDir""$qcStage"multiqc/
