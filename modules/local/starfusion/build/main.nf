process STARFUSION_BUILD {
    tag "$meta.id"
    label 'process_high'

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
    if (dfam_species != "human" && dfam_species != "mouse") {
        error "Invalid species for --dfam_db. Only 'human' or 'mouse' are accepted. Provided: ${dfam_species}"
    }
    def args = task.ext.args ?: ''
    """
    prep_genome_lib.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --dfam_db ${dfam_species} \\
        --pfam_db current \\
        --fusion_annot_lib $fusion_annot_lib \\
        ${args} \\
        --CPU $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ctat_genome_lib_build_dir
    touch ctat_genome_lib_build_dir/AnnotFilterRule.pm
    gzip -c /dev/null > ctat_genome_lib_build_dir/PFAM.domtblout.dat.gz
    touch ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans
    touch ctat_genome_lib_build_dir/ref_genome.fa.mm2
    touch ctat_genome_lib_build_dir/ref_genome.fa.ntf
    gzip -c /dev/null > ctat_genome_lib_build_dir/blast_pairs.dat.gz
    touch ctat_genome_lib_build_dir/ref_annot.cdna.fa
    touch ctat_genome_lib_build_dir/ref_annot.gtf.mini.sortu
    touch ctat_genome_lib_build_dir/ref_genome.fa.ndb
    touch ctat_genome_lib_build_dir/ref_genome.fa.nto
    touch ctat_genome_lib_build_dir/blast_pairs.idx
    touch ctat_genome_lib_build_dir/ref_annot.cdna.fa.idx
    touch ctat_genome_lib_build_dir/ref_annot.gtf.mm2.splice.bed
    touch ctat_genome_lib_build_dir/ref_genome.fa.nhr
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx
    touch ctat_genome_lib_build_dir/__chkpts
    touch ctat_genome_lib_build_dir/ref_annot.cds
    touch ctat_genome_lib_build_dir/ref_annot.pep
    touch ctat_genome_lib_build_dir/ref_genome.fa.nin
    touch ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat
    gzip -c /dev/null > ctat_genome_lib_build_dir/fusion_annot_lib.gz
    touch ctat_genome_lib_build_dir/ref_annot.cdsplus.fa
    touch ctat_genome_lib_build_dir/ref_annot.prot_info.dbm
    touch ctat_genome_lib_build_dir/ref_genome.fa.njs
    touch ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm
    touch ctat_genome_lib_build_dir/fusion_annot_lib.idx
    touch ctat_genome_lib_build_dir/ref_annot.cdsplus.fa.idx
    touch ctat_genome_lib_build_dir/ref_genome.fa
    touch ctat_genome_lib_build_dir/ref_genome.fa.not
    gzip -c /dev/null > ctat_genome_lib_build_dir/trans.blast.dat.gz
    touch ctat_genome_lib_build_dir/pfam_domains.dbm
    touch ctat_genome_lib_build_dir/ref_annot.gtf
    touch ctat_genome_lib_build_dir/ref_genome.fa.fai
    touch ctat_genome_lib_build_dir/ref_genome.fa.nsq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

}
