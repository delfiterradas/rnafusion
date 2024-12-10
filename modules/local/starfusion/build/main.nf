process STARFUSION_BUILD {
    tag 'star-fusion'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/dfam_hmmer_minimap2_samtools_pruned:63e3d21ca68ea531'

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    path fusion_annot_lib
    val dfam_species

    output:
    path "ctat_genome_lib_build_dir"  , emit: reference

    script:
    def binPath = (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1)  ? "prep_genome_lib.pl" : "/opt/conda/lib/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl"
    if (dfam_species != "human" && dfam_species != "mouse") {
        error "Invalid species for --dfam_db. Only 'human' or 'mouse' are accepted. Provided: ${dfam_species}"
    }
    """
    prep_genome_lib.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --dfam_db ${dfam_species} \\
        --fusion_annot_lib $fusion_annot_lib \\
        --CPU $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
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
