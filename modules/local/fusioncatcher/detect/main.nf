process FUSIONCATCHER_DETECT {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/clinicalgenomics/fusioncatcher:1.33"

    input:
    tuple val(meta), path(fastqs, stageAs: "input/*")
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.fusioncatcher.fusion-genes.txt"), emit: fusions, optional: true
    tuple val(meta), path("*.fusioncatcher.summary.txt")     , emit: summary, optional: true
    tuple val(meta), path("*.fusioncatcher.log")             , emit: log
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def single_end = meta.single_end ? "--single-end" : ""
    """
    fusioncatcher.py \\
        -d ${reference} \\
        -i input \\
        -p ${task.cpus} \\
        -o . \\
        --skip-blat \\
        ${single_end} \\
        ${args}

    mv final-list_candidate-fusion-genes.txt ${prefix}.fusioncatcher.fusion-genes.txt
    mv summary_candidate_fusions.txt ${prefix}.fusioncatcher.summary.txt
    mv fusioncatcher.log ${prefix}.fusioncatcher.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher.py --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.fusioncatcher.fusion-genes.txt
    touch ${prefix}.fusioncatcher.summary.txt
    touch ${prefix}.fusioncatcher.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusioncatcher: "\$(fusioncatcher.py --version 2>&1 | awk '{print \$2}')"
    END_VERSIONS
    """
}
