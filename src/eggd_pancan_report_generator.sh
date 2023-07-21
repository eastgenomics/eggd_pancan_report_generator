set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs
    # Move all files in subdirectories of /in directory to the current project
    find ~/in -type f -name "*" -print0 | xargs -0 -I {} mv {} ./

    # Set up samtools to remove chr and create bam

    # samtools in htslib doesn't work as its missing a library, so
    # will install the missing libraries from the downloaded deb files
    # (so as to not use internet)
    sudo dpkg -i libtinfo5_6.2-0ubuntu2_amd64.deb
    sudo dpkg -i libncurses5_6.2-0ubuntu2_amd64.deb

    samtools view -H $bam_name |\
    sed -e 's/chr1/1/' | sed -e 's/chr2/2/' | \
    sed -e 's/chr3/3/' | sed -e 's/chr4/4/' | \
    sed -e 's/chr5/5/' | sed -e 's/chr6/6/' | \
    sed -e 's/chr7/7/' | sed -e 's/chr8/8/' | \
    sed -e 's/chr9/9/' | sed -e 's/chr10/10/' | \
    sed -e 's/chr11/11/' | sed -e 's/chr12/12/' | \
    sed -e 's/chr13/13/' | sed -e 's/chr14/14/' | \
    sed -e 's/chr15/15/' | sed -e 's/chr16/16/' | \
    sed -e 's/chr17/17/' | sed -e 's/chr18/18/' | \
    sed -e 's/chr19/19/' | sed -e 's/chr20/20/' | \
    sed -e 's/chr21/21/' | sed -e 's/chr22/22/' | \
    sed -e 's/chrX/X/' | sed -e 's/chrY/Y/' | \
    sed -e 's/chrM/M/' | samtools reheader - $bam_name > ${bam_prefix}.nochr.bam

    samtools index -b ${bam_prefix}.nochr.bam

    docker load -i $chimerviz_docker_name
    # Get image id from docker image loaded
    CHI_IMAGE_ID=$(sudo docker images --format="{{.Repository}} {{.ID}}" | grep "^sample_report" | cut -d' ' -f2)

    # get the Rscript with the paths in
    echo "Rscript -e \"rmarkdown::render('${rmarkdown_name}', \
    params = list(qcmetrics = '$rnaseqc_metrics_name', \
    qcgenecov = '$rnaseqc_coverage_name', \
    capturegenes = '$capture_bed_name', \
    qcexonscov = '$rnaseqc_exon_name', \
    fusionsabridgedcoding = '$fusioninspector_abridged_name', \
    ensbsqlite = '$ensdb_sqlite_name', \
    bam = '${bam_prefix}.nochr.bam', \
    samplename = 'NA', \
    sexphenotype = 'NA', \
    targetlist = 'NA'))\"" > cmd.sh


    docker run -v /home/dnanexus:/home/software $CHI_IMAGE_ID /bin/bash -c 'ls; cat cmd.sh; bash cmd.sh'

    mkdir -p /home/dnanexus/out/html_report/

    mv *.html /home/dnanexus/out/html_report/
    dx-upload-all-outputs
}
