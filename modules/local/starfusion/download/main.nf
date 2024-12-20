process STARFUSION_DOWNLOAD {
    tag 'star-fusion'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/bed86145102fdf7e381e1a506a4723676f98b4bbe1db5085d02213cef18525c9/data' :
        'community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:aa3a8e3951498552'}"

    output:
    path "ctat_genome_lib_build_dir/*"            , emit: reference
    path "ctat_genome_lib_build_dir/ref_annot.gtf", emit: chrgtf


    script:
    """
    wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz --no-check-certificate

    tar xvf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz

    rm GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz

    mv */ctat_genome_lib_build_dir .
    """

    stub:
    """
    mkdir ctat_genome_lib_build_dir
    touch ref_annot.cdna.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}
