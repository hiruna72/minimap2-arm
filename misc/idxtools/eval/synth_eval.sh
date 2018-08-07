#!/bin/bash


GENERATE_TRUTHSETS(){

MINIMAP=../minimap2-orig/minimap2
PAFTOOLS=../minimap2-orig/misc/paftools.js

#generate indices
echo "Generating index"
$MINIMAP -a -x map-pb ../hg38noAlt.fa -d hg38noAltpb_single.idx > /dev/null 2> hg38noAltpb_single.idx.log
for i in 1700M 850M 350M 140M
do
	echo "Generating index for $i"
    $MINIMAP -a -x map-pb ../hg38noAlt.fa -d hg38noAltpb_$i.idx -I$i > /dev/null 2> hg38noAltpb_$i.idx.log
done
#1700M-2parts 850M-4parts 350M-8parts 140M-16parts

#truthset
echo "Generating SAM truthset"
$MINIMAP -a -x map-pb hg38noAltpb_single.idx pb-1.fa -t 8 > groundtruth.sam 2> groundtruth.sam.log
k8 $PAFTOOLS mapeval groundtruth.sam > groundtruth.sam.eval
echo "Generating PAF truthset"
$MINIMAP -x map-pb hg38noAltpb_single.idx pb-1.fa -t 8 > groundtruth.paf 2> groundtruth.paf.log
k8 $PAFTOOLS mapeval groundtruth.paf > groundtruth.paf.eval

#split index results before modifications to minimap
for i in 1700M 850M 350M 140M
do
	echo "Multi-part index with original minimap2 for $i : SAM"
    $MINIMAP -a -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 > before.$i.sam 2> before.$i.sam.log
    grep -v "^@" before.$i.sam | sort -k1,1  > sorted.before.$i.sam
    k8 $PAFTOOLS mapeval sorted.before.$i.sam >  before.$i.sam.eval
	echo "Multi-part index with original minimap2 for $i : PAF"
    $MINIMAP -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 > before.$i.paf 2> before.$i.paf.log
    sort -k1,1 before.$i.paf > sorted.before.$i.paf
    k8 $PAFTOOLS mapeval sorted.before.$i.paf >  before.$i.paf.eval
done
}

SANITY_CHECK(){
#################After
MINIMAP=../minimap2/minimap2
PAFTOOLS=../minimap2-orig/misc/paftools.js

echo "Generating SAM for whole index"
$MINIMAP -a -x map-pb hg38noAltpb_single.idx pb-1.fa -t 8  > groundtruth_after.sam 2> groundtruth_after.sam.log
diff groundtruth.sam groundtruth_after.sam -q || echo "Diff between before.$i.sam and after.$i.sam failed"
echo "Generating PAF for whole index"
$MINIMAP -x map-pb hg38noAltpb_single.idx pb-1.fa  -t 8 > groundtruth_after.paf 2> groundtruth_after.paf.log
diff groundtruth.paf groundtruth_after.paf -q || echo "Diff between before.$i.sam and after.$i.sam failed"

#split index results after modifications to minimap
for i in 1700M 850M 350M 140M
do
	echo "Multi-part index with modified minimap2 for $i without merging: SAM"
    $MINIMAP -a -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 > after.$i.sam 2> after.$i.sam.log
    diff before.$i.sam after.$i.sam -q || echo "Diff between before.$i.sam and after.$i.sam failed"
	echo "Multi-part index with modified minimap2 for $i without merging: PAF"
    $MINIMAP -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 > after.$i.paf 2> after.$i.paf.log
    diff before.$i.paf after.$i.paf -q || echo "Diff between before.$i.sam and after.$i.sam failed"
done

}

ACCURACY_CHECK(){

MINIMAP=../minimap2/minimap2
PAFTOOLS=../minimap2-orig/misc/paftools.js

#after merge
#split index results after modifications to minimap
for i in 1700M 850M 350M 140M
do
	echo "Multi-part index with modified minimap2 for $i with merging: SAM"
    $MINIMAP -a -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 --multi-prefix tmp > merge.$i.sam 2> merge.$i.sam.log
    k8 $PAFTOOLS mapeval merge.$i.sam >  merge.$i.sam.eval
	echo "Multi-part index with modified minimap2 for $i without merging: PAF"
    $MINIMAP -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 --multi-prefix tmp > merge.$i.paf 2> merge.$i.paf.log
    k8 $PAFTOOLS mapeval merge.$i.paf >  merge.$i.paf.eval
done
}


GENERATE_TRUTHSETS
echo "Truth sets generated"
SANITY_CHECK
echo "Sanity Check completed"
ACCURACY_CHECK
echo "Accuracy check completed"
