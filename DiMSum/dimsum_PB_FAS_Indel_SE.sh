#!/bin/bash

#Run DiMSum Stages 1-5 in single-end mode using only demultiplexed forward reads that have full coverage of the exon sequence

fastqFileDir="PB_FAS_Indel/tmp/0_demultiplex"
gzipped="FALSE"
experimentDesignPath="experimentDesign_PB_FAS_Indel_SE.txt"
cutadaptMinLength="3"
cutadaptErrorRate="0.2"
vsearchMinQual="30"
vsearchMaxee="0.5"
outputPath="./"
projectName="PB_FAS_Indel_SE"
wildtypeSequence="GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG"
numCores="10"
maxSubstitutions="63"
experimentDesignPairDuplicates="T"
cutadaptOverlap="1"
indels="all"
retainIntermediateFiles="T"
paired="F"
barcodeIdentityPath="variantIdentity_PB_FAS_Indel.txt"

DiMSum -i $fastqFileDir -g $gzipped -e $experimentDesignPath -n $cutadaptMinLength -q $vsearchMinQual -m $vsearchMaxee -o $outputPath -p $projectName -w $wildtypeSequence -c $numCores --maxSubstitutions $maxSubstitutions --experimentDesignPairDuplicates $experimentDesignPairDuplicates --cutadaptOverlap $cutadaptOverlap --indels $indels --retainIntermediateFiles $retainIntermediateFiles --paired $paired --barcodeIdentityPath $barcodeIdentityPath
