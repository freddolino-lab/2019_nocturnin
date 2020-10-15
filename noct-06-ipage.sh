#!/usr/bin/env bash
set -e

ANALYSISPATH="$HOME"'/data/noct/humanOE/revision/z-iPAGE/'
# # run 1: iPAGE v1.2a from /home/diaorch/src/iPAGE
# INPUTFILENAME="$ANALYSISPATH"'20200604-noct-z_shrunken_input.tsv'
# PAGEDBNAME='noct201806'
# perl ${PAGEDIR}/page.pl --expfile="$INPUTFILENAME" --species="$PAGEDBNAME" --exptype=continuous --ebins 9 > "$ANALYSISPATH"'z_shrunken_ipage.bin9.log' 2> "$ANALYSISPATH"'z_shrunken_ipage.bin9.err'

# # run 2: iPAGE v1.2a from /home/petefred/src/iPAGE
# # the two input files are identical (the 20200605 file is a symbolic link to the 20200604 one), they were run separately to use the iPAGE version from /home/petefred/src/iPAGE so the expression-bin-significance bold box would show up
# # according to iPAGE/page.pl, they usage is:
# # "Usage: perl page.pl --expfile=FILE --datafile=FILE --goindexfile=FILE --gonamesfile=FILE --exptype=TXT --cattypes=F,C,P --catmaxcount=INT\n"
# # combine with exmaple in ~/data/noct/humanOE/revision/z-iPAGE/20200604-noct-z_shrunken_input.tsv_PAGE/info.txt
INPUTFILENAME="$ANALYSISPATH"'20200605-noct-z_shrunken_input.tsv.sigBox'
myPAGEDIR="$PAGEDIR"
PAGEDIR='/home/petefred/src/iPAGE'
GODB='/home/diaorch/src/iPAGE/PAGE_DATA/ANNOTATIONS/noct201806/noct201806'
# perl ${PAGEDIR}/page.pl --expfile="$INPUTFILENAME" --goindexfile="$GODB"'_index.txt' --gonamesfile="$GODB"'_names.txt' --exptype=continuous --ebins 9 > "$ANALYSISPATH"'z_shrunken_ipage.bin9.sigBox.log' 2> "$ANALYSISPATH"'z_shrunken_ipage.bin9.sigBox.err'
# PAGEDIR="$myPAGEDIR"
INPUT_INDP0="$INPUTFILENAME"'.indp0'
ln -s "$INPUTFILENAME" "$INPUT_INDP0"
perl ${PAGEDIR}/page.pl --expfile="$INPUT_INDP0" --goindexfile="$GODB"'_index.txt' --gonamesfile="$GODB"'_names.txt' --exptype=continuous --ebins 9 --independence=0 > "$ANALYSISPATH"'z_shrunken_ipage.bin9.sigBox.indp0.log' 2> "$ANALYSISPATH"'z_shrunken_ipage.bin9.sigBox.indp0.err'
