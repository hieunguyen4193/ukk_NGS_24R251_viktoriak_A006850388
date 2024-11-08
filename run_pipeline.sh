main_nf="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/rnaseq/main.nf";
SampleSheet="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/SampleSheet.csv";
outdir="/media/hieunguyen/HNSD01/outdir/NGS_24R251_viktoriak_A006850388";
mkdir -p $outdir;

work="./work";

##### mouse reference
# fasta="/home/hieunguyen/CRC1382/src_2023/RNAseq_aT/annotation_resources/Mus_musculus.GRCm39.dna.primary_assembly.fa";
# gtf="/home/hieunguyen/CRC1382/src_2023/RNAseq_aT/annotation_resources/Mus_musculus.GRCm39.109.gtf";

##### human reference
# fasta="/home/hieunguyen/CRC1382/storage/reference_genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
# gtf="/home/hieunguyen/CRC1382/storage/reference_genomes/Homo_sapiens.GRCh38.108.gtf.gz";

nextflow="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/nextflow";

max_mem="75GB";
max_cpus=20;

${nextflow} run ${main_nf} --input $SampleSheet --outdir ${outdir} --genome GRCh37 -profile docker -resume --max_cpus $max_cpus  --max_memory $max_mem --publish_dir_mode copy
