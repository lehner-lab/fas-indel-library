#!/bin/bash

#Run DiMSum in default paired-end mode to demultiplex reads into input and replicate output samples (Stage 0 only)

fastqFileDir="PATH_TO_FASTQ_DIR"
experimentDesignPath="experimentDesign_PB_FAS_Indel.txt"
barcodeDesignPath="barcodeDesign_PB_FAS_Indel.txt"
outputPath="./"
projectName="PB_FAS_Indel"
wildtypeSequence="GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG"
stranded="F"
experimentDesignPairDuplicates="T"
startStage="0"
stopStage="0"
retainIntermediateFiles="T"

DiMSum -i $fastqFileDir -e $experimentDesignPath -b $barcodeDesignPath -o $outputPath -p $projectName -w $wildtypeSequence --stranded $stranded --experimentDesignPairDuplicates $experimentDesignPairDuplicates -s $startStage --stopStage $stopStage --retainIntermediateFiles $retainIntermediateFiles
