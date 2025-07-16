process FUSIONCATCHER_DOWNLOAD {
    tag "fusioncatcher_download"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b97801b2147f99c974bdf01b821c043654f6bf1f2864a1633231740999df072/data':
        'community.wave.seqera.io/library/tar_pip_gdown:34715940085a0f1c' }"

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
    gdown $args --fuzzy ${params.fusioncatcher_download_link}
    tar xz human_v${genome_gencode_version}.tar.gz
    rm human_v${genome_gencode_version}.tar*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gdown: \$(gdown --version | head -1 | cut -d ' ' -f 2)
        tar: \$(tar --version | head -1 | sed -e 's/tar (GNU tar) //')
    END_VERSIONS
    """

    stub:
    """
    mkdir human_v${genome_gencode_version}
    touch human_v${genome_gencode_version}/ensembl_fully_overlapping_genes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gdown: \$(gdown --version | head -1 | cut -d ' ' -f 2)
        tar: \$(tar --version | head -1 | sed -e 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
