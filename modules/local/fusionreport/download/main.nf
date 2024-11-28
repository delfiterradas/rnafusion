process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    conda "bioconda::star=2.7.9a"
    container "docker.io/clinicalgenomics/fusion-report:3.1.0"

    input:
    val(username)
    val(passwd)

    output:
    path "*"                , emit: reference
    path "versions.yml"     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: ''
    """
    fusion_report download --cosmic_usr "$username" --cosmic_passwd "$passwd" $args ./

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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """

}
