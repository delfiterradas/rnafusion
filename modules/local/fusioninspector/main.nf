process FUSIONINSPECTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d98081910b41caa0b62d892e7a44543183111afa2b0a5eea5a4a9472eea72611/data' :
        'community.wave.seqera.io/library/fusion-inspector:2.10.0--b59230553c98db81'}"

    input:
    tuple val(meta), path(reads), path(fusion_list)
    path reference

    output:
    tuple val(meta), path("*FusionInspector.fusions.tsv")                  , emit: tsv
    tuple val(meta), path("*.coding_effect")                , optional:true, emit: tsv_coding_effect
    tuple val(meta), path("*.gtf")                          , optional:true, emit: out_gtf
    path "*"                                                               , emit: output
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    FusionInspector \\
        --fusions $fusion_list \\
        --genome_lib ${reference} \\
        $fasta \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix $prefix \\
        --vis $args $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.FusionInspector.fusions.tsv.annotated.coding_effect
    touch ${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}
