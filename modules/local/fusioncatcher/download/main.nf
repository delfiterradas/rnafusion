process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fusioncatcher:1.33--hdfd78af_5':
        'biocontainers/fusioncatcher:1.33--hdfd78af_5' }"

    input:
    val genome_gencode_version

    output:
    tuple env(meta), path("human_v*"), emit: reference
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    meta = [ id: "human_v${genome_gencode_version}" ]
    """
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partaa
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partab
    wget --no-check-certificate $args ${params.fusioncatcher_download_link}/human_v${genome_gencode_version}.tar.gz.partac
    cat human_v${genome_gencode_version}.tar.gz.* | tar xz
    rm human_v${genome_gencode_version}.tar*

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
