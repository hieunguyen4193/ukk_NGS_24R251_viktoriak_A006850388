main_nf="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/rnaseq/main.nf";
# SampleSheet="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/SampleSheet.csv";
# SampleSheet="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/LYN_CSK_samplesheet.csv"
SampleSheet="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/LYN_CSK_samplesheet_bosung.csv"

outdir="/media/hieunguyen/HNSD01/outdir/NGS_24R251_viktoriak_A006850388_replace";
work="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/work"

mkdir -p ${work}
mkdir -p $outdir;

##### mouse reference
# fasta="/home/hieunguyen/CRC1382/src_2023/RNAseq_aT/annotation_resources/Mus_musculus.GRCm39.dna.primary_assembly.fa";
# gtf="/home/hieunguyen/CRC1382/src_2023/RNAseq_aT/annotation_resources/Mus_musculus.GRCm39.109.gtf";

##### human reference
# fasta="/home/hieunguyen/CRC1382/storage/reference_genomes/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
# gtf="/home/hieunguyen/CRC1382/storage/reference_genomes/Homo_sapiens.GRCh38.108.gtf.gz";

nextflow="/media/hieunguyen/HNSD01/src/NGS_24R251_viktoriak_A006850388/nextflow";
path_to_ref="/media/hieunguyen/HNHD01/storage/ref";

# gtf=${path_to_ref}/Mus_musculus.GRCm39.113.gtf;
# fasta=${path_to_ref}/Mus_musculus.GRCm39.dna.primary_assembly.fa;

gtf=${path_to_ref}/Homo_sapiens.GRCh38.108.gtf;
fasta=${path_to_ref}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa;
max_mem="75GB";
max_cpus=20;


${nextflow} run ${main_nf} \
--input $SampleSheet \
--outdir ${outdir} \
--fasta $fasta --gtf $gtf \
-profile docker \
-resume --max_cpus $max_cpus  --max_memory $max_mem --publish_dir_mode copy -w ${work}
