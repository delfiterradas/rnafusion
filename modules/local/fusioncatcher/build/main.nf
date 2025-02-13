process FUSIONCATCHER_BUILD {
    tag "fusioncatcher_build"
    label 'process_medium'

    container "docker.io/clinicalgenomics/fusioncatcher:1.33"

    input:
    val genome_gencode_version

    output:
    tuple env(meta), path("human_v${genome_gencode_version}"), emit: reference
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    meta = [ id: "human_v${genome_gencode_version}" ]
    """
    fusioncatcher-build.py \\
        -g homo_sapiens \\
        -o human_v${genome_gencode_version} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
