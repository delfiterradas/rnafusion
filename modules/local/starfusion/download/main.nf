process STARFUSION_DOWNLOAD {
    tag 'star-fusion'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/dfam_hmmer_samtools_star-fusion_pruned:5694d82381bf039e'

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
