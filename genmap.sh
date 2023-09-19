conda install -c bioconda genmap
conda activate genmap
cd /groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2
genmap index -F Anthoceros_agrestis_Oxford_genome.fasta -I genmapid
genmap map -K 50 -E 0 -I genmapidx/ -O /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap -w -bg
cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq
sort -k1,1 -k2,2n genmap.bedgraph > sorted.genmap.bedgraph

module load build-env/.f2021
module load  bedgraphtobigwig/385-linux-x86_64
bedGraphToBigWig sorted.genmap.bedgraph /groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Anthoceros_agrestis_Oxford_chrom_sizes.txt AagrOXFv2_genmap_k50_E0.bw
