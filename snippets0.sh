PROJECT_DIR=/Volumes/silver/phd/20231217/ice_water


declare -a NAMES=(
    'Water-2a'
    'Water-3a'
    'Water-5a'
    'Water-6a'
    'Water-7a'
    'Water-9a'
    'Water-10a'
    'Water-Neg1'
    'Water-Neg2'
    'Water-Neg3'
)

# # for fastp trimming
# cd ${PROJECT_DIR}
# for name in "${NAMES[@]}" ; do
#     echo "$name"
#     fastp -i raw/${name}_1.fastq.gz -I raw/${name}_2.fastq.gz -o trimmed/${name}_1.trimmed.fq -O trimmed/${name}_2.trimmed.fq --json fastp/${name}.json --html fastp/${name}.html > fastp/${name}.log
# done

############  QIIME  ############

micromamba activate qiime2-2023.7
PROJECT_DIR=/Volumes/silver/phd/20231217/ice_water
DATA_DIR="${PROJECT_DIR}/trimmed"
SUFFIX=.trimmed
FQ_SUFFIX=.trimmed.fq

# to generate qiime2 manifest
MANIFEST=$PROJECT_DIR/manifest${SUFFIX}.tsv

printf '%s\t%s\t%s\n' \
    'sample-id' \
    'forward-absolute-filepath' \
    'reverse-absolute-filepath' \
    > "${MANIFEST}"

for name in "${NAMES[@]}" ; do
    printf '%s\t%s\t%s\n' \
        "$name" \
        "${DATA_DIR}/${name}_1.${FQ_SUFFIX}" \
        "${DATA_DIR}/${name}_2.${FQ_SUFFIX}" \
        >> "${MANIFEST}"
done

PROJECT_DIR=/Volumes/silver/phd/20231217/ice_water
SUFFIX=.trimmed
MANIFEST=$PROJECT_DIR/manifest${SUFFIX}.tsv
PE_QZA=demux-PE${SUFFIX}.qza

###  1. Import qiime tools. GZIPPED FILES ARE MUCH FASTER  ##########
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path "${MANIFEST}" \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path ${PE_QZA}




##########  2a. Denoise, dereplicate and filters chimeras with DADA2  ##########
# SUFFIX=.trimmed
SUFFIX=
TRIM_LEN=30
TRUNC_LEN=0
THREADS=4
# input
DEMUX=demux-PE${SUFFIX}.qza
# output
DADA2_TABLE="dada2_table${SUFFIX}.qza"
DADA2_REP_SEQS="dada2_rep_seqs${SUFFIX}.qza"
DADA2_DENOISE_STATS="dada2_denoise_stats${SUFFIX}.qza"

time qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "${DEMUX}" \
    --p-trim-left-f ${TRIM_LEN} \
    --p-trim-left-r ${TRIM_LEN} \
    --p-trunc-len-f ${TRUNC_LEN} \
    --p-trunc-len-r ${TRUNC_LEN} \
    --p-n-threads ${THREADS} \
    --verbose \
    --o-table "${DADA2_REP_TABLE}" \
    --o-representative-sequences "${DADA2_REP_SEQS}" \
    --o-denoising-stats "${DADA2_DENOISE_STATS}" \
    1>dada2_denoise_paired${SUFFIX}.log 2>&1



## Feature table summary
# Viz: make metadata file in tsv format: i.e. barcode_metadata.tsv
## Generate visual and tabular summaries of a feature table.
PROJECT_DIR=/Users/vivianleung/Dropbox/phd/ice_water
SUFFIX=
# input
INPUT_TABLE="dada2_table${SUFFIX}.qza"
METADATA="${PROJECT_DIR}/metadata.tsv"
# output
VIZ_FILE="${INPUT_TABLE/%a/v}"

time qiime feature-table summarize \
    --i-table "${INPUT_TABLE}" \
    --o-visualization "${VIZ_FILE}" \
    --m-sample-metadata-file "${METADATA}"


##########  2. Quality filter by Q-score  ##########
SUFFIX= # .trimmed
# input
DEMUX=demux-merged${SUFFIX}.qza
# output
FILTERED_SEQS="filtered${SUFFIX}.qza"
FILTER_STATS="filter_stats${SUFFIX}.qza"

time qiime quality-filter q-score \
--i-demux "${DEMUX}" \
--o-filtered-sequences "${FILTERED_SEQS}" \
--o-filter-stats "${FILTER_STATS}"

# convert back to PairedEndSequencesWithQuality for merge-pairs
SUFFIX= # .trimmed
INPUT_PATH=demux-filtered${SUFFIX}.qza
OUTPUT_PATH=demux-filtered-PE${SUFFIX}.qza
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "${INPUT_PATH}" \
    --output-path "${OUTPUT_PATH}" 


## Viz: Visualize filter results
# SUFFIX=.trimmed
SUFFIX=
INPUT_FILE="filter_stats${SUFFIX}.qza"
VIZ_FILE="${INPUT_FILE/%a/v}"

time qiime metadata tabulate \
 --m-input-file "${INPUT_FILE}" \
 --o-visualization "${VIZ_FILTER_STATS}"
 

# Viz: Generate tabular view of feature identifier to sequence mapping, 
#   including links to BLAST each sequence against the NCBI nt database.
SUFFIX=
INPUT_SEQS="dada2_rep_seqs${SUFFIX}.qza"
VIZ_FILE="${INPUT_SEQS/%a/v}"

qiime feature-table tabulate-seqs \
    --i-data ${INPUT_SEQS} \
    --o-visualization ${VIZ_FILE}


# Viz: Generate a tabular view of Metadata. The output visualization supports 
#   interactive filtering, sorting, and exporting to common file formats.
SUFFIX=
INPUT_FILE="dada2_denoise_stats${SUFFIX}.qza"
VIZ_FILE="${INPUT_FILE/%a/v}"

qiime metadata tabulate \
    --m-input-file "${INPUT_FILE}" \
    --o-visualization "${VIZ_FILE}"


##########  3b.1 Merge paired end reads with vsearch  ##########
# SUFFIX=.trimmed
SUFFIX=
DEMUX=demux${SUFFIX}.qza
MERGED_DEMUX=demux-merged${SUFFIX}.qza
THREADS=8

date && time qiime vsearch merge-pairs \
    --p-threads $THREADS \
    --p-minlen 50 \
    --verbose \
    --i-demultiplexed-seqs "${DEMUX}" \
    --o-merged-sequences "${MERGED_DEMUX}"


## Viz: check merge results for trim and trunc lengths to apply
SUFFIX=
INPUT_DATA=filtered-merged${SUFFIX}.qza
VIZ_FILE="${INPUT_DATA/%a/v}"
qiime demux summarize \
    --i-data "${INPUT_DATA}" \
    --o-visualization "${VIZ_FILE}"

########## 3b.2 Seq control with Deblur on filtered data  ##########
# SUFFIX=.trimmed
SUFFIX=
TRIM_LEN=30
LEFT_TRIM_LEN=30
# input
DEMUX="filtered-merged${SUFFIX}.qza"
# output
DEBLUR_REP_SEQS="deblur_rep_seqs${SUFFIX}.qza"
DEBLUR_TABLE="deblur_table${SUFFIX}.qza"
DEBLUR_STATS="deblur_stats${SUFFIX}.qza"

qiime deblur denoise-16S \
 --p-trim-length ${TRIM_LENGTH}  \
 --p-left-trim-len ${LEFT_TRIM_LEN}  \
 --i-demultiplexed-seqs "${DEMUX}" \
 --o-representative-sequences "${DEBLUR_REP_SEQS}" \
 --o-table "${DEBLUR_TABLE}" \
 --p-sample-stats \
 --o-stats "${DEBLUR_STATS}"


# Viz: visualize deblur stats
SUFFIX=.raw
INPUT_DEBLUR_STATS="deblur_stats${SUFFIX}.qza"
VIZ_FILE="${INPUT_DEBLUR_STATS/%a/v}"
qiime deblur visualize-stats \
 --i-deblur-stats "${INPUT_DEBLUR_STATS}" \
 --o-visualization "${VIZ_FILE}"


##########  Phylogeny  ##########
SUFFIX=
# input 
DEBLUR_REP_SEQS="deblur_rep_seqs${SUFFIX}.qza"
# output
ALIGNED_REP_SEQS="deblur.aligned${SUFFIX}.qza"
MASKED_ALIGNED_REP_SEQS="deblur.masked_aligned${SUFFIX}.qza"
UNROOTED_TREE="deblur.unrooted_tree${SUFFIX}.qza"
ROOTED_TREE="deblur.rooted_tree${SUFFIX}.qza"
THREADS=4
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "${DEBLUR_REP_SEQS}" \
    --p-n-threads ${THREADS} \
    --verbose \
    --o-alignment "${ALIGNED_REP_SEQS}" \
    --o-masked-alignment "${MASKED_ALIGNED_REP_SEQS}" \
    --o-tree "${UNROOTED_TREE}" \
    --o-rooted-tree "${ROOTED_TREE}"


##########  Classification  ##########
PROJECT_DIR=/Users/vivianleung/Dropbox/phd/ice_water
SUFFIX=
THREADS=4
# input
DEBLUR_REP_SEQS="deblur_rep_seqs${SUFFIX}.qza"
REF_READS=${PROJECT_DIR}/external/silva-138-99-seqs.qza
REF_TAXONOMY=${PROJECT_DIR}/external/silva-138-99-tax.qza
# output
CLASSIFICATION="deblur.taxonomy${SUFFIX}.qza"

qiime feature-classifier classify-consensus-vsearch \
 --i-query "${DEBLUR_REP_SEQS}" \
 --i-reference-reads "${REF_READS}" \
 --i-reference-taxonomy "${REF_TAXONOMY}" \
 --o-classification "${CLASSIFICATION}" \
 --p-perc-identity 0.97 \
 --p-maxaccepts 3 \
 --p-threads $THREADS

## Viz: Visualize classification results
SUFFIX=
INPUT_FILE="taxonomy${SUFFIX}.qza"
VIZ_FILE="${INPUT_FILE/%a/v}"

qiime metadata tabulate \
--m-input-file ${INPUT_FILE} \
--o-visualization ${VIZ_FILE}

## Viz: Visualize classification results
SUFFIX=
# input
DEBLUR_TABLE="deblur_table${SUFFIX}.qza"
INPUT_CLASSIFICATION="deblur.taxonomy${SUFFIX}.qza"
METADATA="${PROJECT_DIR}/metadata.tsv"
# output
VIZ_FILE="${INPUT_CLASSIFICATION/%qza/boxplot.qzv}"

qiime taxa barplot \
    --i-table "${DEBLUR_TABLE}" \
    --i-taxonomy "${INPUT_CLASSIFICATION}" \
    --m-metadata-file "${METADATA}" \
    --o-visualization "${VIZ_FILE}"



##########  Diversity analysis  ##########
SUFFIX=
SAMPLING_DEPTH=80000
# input
ROOTED_TREE="deblur.rooted_tree${SUFFIX}.qza"
DEBLUR_TABLE="deblur_table${SUFFIX}.qza"
PROJECT_DIR=/Users/vivianleung/Dropbox/phd/ice_water
METADATA="${PROJECT_DIR}/metadata.tsv"

# output
DIVERSITY_MERGED_RESULTS="diversity_metrics_results"

# mkdir -p "${DIVERSITY_MERGED_RESULTS}"

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "${ROOTED_TREE}" \
    --i-table "${DEBLUR_TABLE}" \
    --m-metadata-file "${METADATA}" \
    --p-sampling-depth ${SAMPLING_DEPTH} \
    --output-dir "${DIVERSITY_MERGED_RESULTS}"

