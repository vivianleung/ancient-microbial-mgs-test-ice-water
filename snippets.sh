#!/usr/bin/env bash

# snippets



######################################################


# To do:
# Kraken2





# Metaphlan
# Kaiju

######################################################

# 30 Dec 2023

###### dada2 of trimmed_PE
qiime tools import --type SampleData[PairedEndSequencesWithQuality] --input-path manifest.trimmed.tsv --input-type PairedEndFastqManifestPhred33V2


qiime dada2 denoise-paired --i-demultiplexed-seqs demux-trimmed-PE.qza --p-trim-left-f 30 --p-trim-left-r 30 --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads 4 --verbose --o-table dada2_denoise-table.qza --o-representative-sequences demux-dada2_denoise-rep_seqs.qza --o-denoising-stats dada2_denoise-stats.qza 1>dada2_denoise.log 2>&1

### >>> RESUME HERE
mv demux-dada2_denoise-rep_seqs.qza dada2_denoise-rep_seqs.qza

qiime feature-table tabulate-seqs \
    --i-data dada2_denoise-rep_seqs.qza \
    --o-visualization dada2_denoise-rep_seqs.qzv

qiime metadata tabulate \
    --m-input-file dada2_denoise-stats.qza \
    --o-visualization dada2_denoise-stats.qzv


# scp -r d24h:/home/d24h_prog2/data_share/ice_water/qiime/{dada2_denoise-rep_seqs.qzv,dada2_denoise-stats.qzv} .

## 

# qiime phylogeny align-to-tree-mafft-fasttree \
#     --i-sequences dada2_denoise-rep_seqs.qza \
#     --p-n-threads 4 \
#     --verbose \
#     --o-alignment aligned.dada2.qza \
#     --o-masked-alignment masked_aligned.dada2.qza \
#     --o-tree unrooted.dada2.qza \
#     --o-rooted-tree rooted.dada2.qza




###### deblur of merged only
qiime tools import --type SampleData[JoinedSequencesWithQuality] \
    --input-path ../merged/manifest.merged.tsv \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path demux-trimmed-merged.qza


qiime tools peek demux-trimmed-merged.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-trimmed-merged.qza \
  --p-trim-length 30 \
  --p-left-trim-len 30 \
  --p-sample-stats \
  --o-representative-sequences demux-deblur_denoise-rep_seqs.merged.qza \
  --o-table deblur_denoise-table.merged.qza \
  --o-stats deblur_denoise-stats.merged.qza


###### deblur of merged_fwd

qiime tools import --type SampleData[SequencesWithQuality]  --input-path ../.
./merged/merged_fwd/manifest.merged_fwd.tsv --input-format SingleEndFastqManifestPhred33V2 --output-path demux-trimmed_merged_fwd.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-trimmed_merged_fwd.qza \
  --p-trim-length 30 \
  --p-left-trim-len 30 \
  --p-sample-stats \
  --o-representative-sequences demux-deblur_denoise-rep_seqs.merged_fwd.qza \
  --o-table deblur_denoise-table.merged_fwd.qza \
  --o-stats deblur_denoise-stats.merged_fwd.qza

### >>> RESUME HERE
mv demux-deblur_denoise-rep_seqs.merged_fwd.qza deblur_denoise-rep_seqs.merged_fwd.qza

qiime deblur visualize-stats \
 --i-deblur-stats deblur_denoise-stats.merged_fwd.qza \
 --o-visualization deblur_denoise-stats.merged_fwd.qzv

qiime feature-table tabulate-seqs \
    --i-data deblur_denoise-rep_seqs.merged_fwd.qza \
    --o-visualization deblur_denoise-rep_seqs.merged_fwd.qzv

# scp -r d24h:/home/d24h_prog2/data_share/ice_water/qiime/merged_fwd/{deblur_denoise-stats.merged_fwd.qzv,deblur_denoise-rep_seqs.merged_fwd.qzv} .

##

# qiime phylogeny align-to-tree-mafft-fasttree \
#     --i-sequences deblur_denoise-rep_seqs.merged_fwd.qza \
#     --p-n-threads 4 \
#     --verbose \
#     --o-alignment aligned.deblur.merged_fwd.qza \
#     --o-masked-alignment masked_aligned.deblur.deblur.merged_fwd.qza \
#     --o-tree unrooted.deblur.merged_fwd.qza \
#     --o-rooted-tree rooted.deblur.merged_fwd.qza




