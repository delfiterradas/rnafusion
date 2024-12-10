process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "bioconda::star=2.7.9a"
    container "docker.io/clinicalgenomics/fusion-report:3.1.0"

    output:
    tuple val(meta), path("fusionreport_dbs"), emit: fusionreport_db
    path "versions.yml"                      , emit: versions

    script:
    meta = [id: 'fusionreport_dbs']
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    fusion_report download $args ./
    mkdir fusionreport_dbs
    mv *.txt *.log *.db fusionreport_dbs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """

    stub:
    meta = [id: 'fusionreport_dbs']
    """
    mkdir fusionreport_dbs
    touch fusionreport_dbs/cosmic.db
    touch fusionreport_dbs/fusiongdb2.db
    touch fusionreport_dbs/mitelman.db
    touch fusionreport_dbs/DB-timestamp.txt
    touch fusionreport_dbs/fusion_report.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
