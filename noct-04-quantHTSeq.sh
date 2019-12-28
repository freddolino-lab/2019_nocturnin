#/usr/bin/env bash

# Re-quantification by HTSeq, specifying stranded as argument

# set-up environment
## Python3 environment with dependencies
## variable dependency
## make sure "$project"/noct-00-envSetup.sh is available  

source noct-00-envSetup.sh

mkdir "$htDir"

# run HTseq
for filename in "$starDir"*"$starSortedBamSuffix" 
do
    # echo "$filename"
    filenameBase=${filename#$starDir}
    filenameBase=${filenameBase%$starSortedBamSuffix}
    echo "$filenameBase"
    "$htPath" --format=bam --order=pos --stranded=reverse --minaqual=10 --mode=union --nonunique=none "$filename" "$newRefName".gtf >"$htDir""$filenameBase""$htSuffix" 2>"$htDir""$filenameBase"'.err'
done

