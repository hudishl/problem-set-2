


#! /usr/bin/env bash


# Use BEDtools intersect to identify the size of the largest overlap
#between CTCF and H3K4me3 locations.

datasets=$HOME/data-sets
bed=$datasets/bed

answer_1=$(bedtools intersect -a $bed/encode.h3k4me3.hela.chr22.bed.gz \
            -b $bed/encode.tfbs.chr22.bed.gz \
    | awk 'BEGIN {OFS="\t"} {print $0, $3-$2}' \
    | sort -k11nr | cut -f11 | head -n1)

echo "answer_1: $answer_1"

echo -e  "chr22\t19000000\t19000500" > interval.bed
answer_2=$(bedtools nuc -fi $datasets/fasta/hg19.chr22.fa -bed $Desktop/GAW/problem-set-2/interval.bed \
        | cut -f5 \
        | tail -n1)

        echo "answer-2: $answer_2"

#Use BEDtools to identify the length of the CTCF ChIP-seq peak
#(i.e., interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

gzcat $bed/encode.tfbs.chr22.bed.gz | grep -w 'CTCF' | gzip -c > ctcfpeaks.bed.gz
answer_3=$(bedtools map -c 4 -o mean -a ctcfpeaks.bed.gz \
            -b $datasets/bedtools/ctcf.hela.chr22.bg.gz \
           | sort -k5nr | head -n1)

echo "answer_3": $answer_3

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of a TSS) with the
#highest median signal in ctcf.hela.chr22.bg.gz. Report the gene name (e.g., 'ABC123')

bedtools slop -b 1000 -i $bed/tss.hg19.chr22.bed.gz \
    -g $datasets/bedtools/hg19.genome > promoters.22.slop1kb.bed
answer_4=$(bedtools map -c 4 -o mean -a promoters.22.slop1kb.bed \
             -b $datasets/bedtools/ctcf.hela.chr22.bg.gz \
            | sort -k7nr | head -n1 | cut -f4)

echo "Answer_4": $answer_4
