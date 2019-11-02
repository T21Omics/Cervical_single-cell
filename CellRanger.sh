set -e
umask 0002
lib = LIB_NAME
sub = SUB_NAME
cellnum=6000
current=$PWD
line1=DATA_PATH
export PATH = /software/cellranger-3.0.2/cellranger-cs/3.0.2/bin:$PATH
cellranger count --id = ${lib}_$sub \
--localcores = 24 \
--localmem = 400 \
--expect-cell = $cellnum \
--transcriptome = GRCh38.p12_index \
--fastqs = $line1 
#remove intermediate files
echo "rm -rf $current/${lib}_$sub/SC_RNA_COUNTER_CS" > $current/rm_SC_RNA_COUNTER_CS.sh
