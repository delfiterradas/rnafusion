process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "bioconda::star=2.7.9a"
    container "docker.io/clinicalgenomics/fusion-report:3.1.0"

    output:
    path "fusiongdb2.db"    , emit: fusiongdb2
    path "mitelman.db"      , emit: mitelman
    path "cosmic.db"        , emit: cosmic, optional: true
    path "*.txt"            , emit: timestamp
    path "*.log"            , emit: log
    path "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''
    """
    fusion_report download $args ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """

    stub:
    """
    touch cosmic.db
    touch fusiongdb2.db
    touch mitelman.db
    touch DB-timestamp.txt
    touch fusion_report.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """

}
