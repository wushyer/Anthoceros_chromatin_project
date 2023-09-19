This is a document describing how we aligned and analyzed ChIP-seq, RNA-seq and EM-seq data presented in

Hisanaga et al., 2023, The ancestral chromatin landscape of land plants, New Phytologist.

===============ChIPseq data analysis==================
ChIP-seq reads were preprocessed and aligned to A. agrestis genome using "ChIPseq-align.sh"

Pearson correlation matrices were generated using "multibamsum_plotcor.sh"

Peaks were called by using "macs2-bamcompare.sh" and overlap between each biological replicate were calculated using "calcoverlap.sh"

Deduplicated reads from 2 biological replicates were merged using "mergeBam"

Using merged bam files above, bigwig files of each mark normalized by H3 were generated using "macs2-bamcompare.sh" 

Peak call was done by using "macs2-bamcompare.sh"

overlap between peaks and genomic features were summarized using "summarize_feature_overlap.sh" (Fig.1a)

Clustering analysis was done by using "computematrix-heatmap.sh"

Optimal numbers of clusters in the k-means clustring was determined by using "Elbow_method2.R"

Cluster names were modified using "modify_cluster_names.R"

Overlaps between TEs and PCGs were calculated using "Calc_TE-PCG_overlap_perCluster.sh"

Genome mappability was calculated and an output bedgraph was converted to a bigwig file using "genmap.sh"

Mappability per each TE or PCG was calculated using "bigwigaverageoverbed_mappability.sh"

closest features of each TE or PCG were identified using "closestfeatures.sh"

plofile plots of H3K4me3 and H3K36me3 over PCGs in cluster P3 and P4 were plotted using "FigS2e_protprofile.sh"

===============RNAseq data analysis==================

RNAseq reads were aligned and count data wwere obtained using "rsem_mapping_table.sh" and "rsem_run.sh"

===============EMseq data analysis===================

EMseq reads were aligned and methlylation data were extracted using "bismark.sh"

methylation information was summarized to bigwig files using "context_to_bw_merged.sh"

DNA methylation levels are plotted over each genomic feature using "plotPlofileDNAmet.sh"


===============plots by R scripts===================

Plots in main figures and supplemental figures were generated using "plots_main.R" and "plots_sup.R", respectively