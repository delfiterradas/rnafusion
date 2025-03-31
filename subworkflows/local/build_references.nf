/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { GENCODE_DOWNLOAD }                from '../../modules/local/gencode_download/main'
include { FUSIONCATCHER_BUILD }             from '../../modules/local/fusioncatcher/build/main'
include { FUSIONREPORT_DOWNLOAD }           from '../../modules/local/fusionreport/download/main'
include { HGNC_DOWNLOAD }                   from '../../modules/local/hgnc/main'
include { STARFUSION_BUILD }                from '../../modules/local/starfusion/build/main'
include { GTF_TO_REFFLAT }                  from '../../modules/local/uscs/custom_gtftogenepred/main'
include { GET_RRNA_TRANSCRIPTS }            from '../../modules/local/get_rrna_transcript/main'
include { CTATSPLICING_PREPGENOMELIB }      from '../../modules/local/ctatsplicing/prepgenomelib/main.nf'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
include { ARRIBA_DOWNLOAD }                 from '../../modules/nf-core/arriba/download/main'
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE }             from '../../modules/nf-core/star/genomegenerate/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST }         from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { SALMON_INDEX }                    from '../../modules/nf-core/salmon/index/main'
include { GFFREAD }                         from '../../modules/nf-core/gffread/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow BUILD_REFERENCES {

    take:
    tools // list of all the tools to create references from

    main:
    ch_versions = Channel.empty()

    if (!exists_not_empty(params.fasta) || !exists_not_empty(params.gtf)){
        GENCODE_DOWNLOAD(params.genome_gencode_version, params.genome)
        ch_versions = ch_versions.mix(GENCODE_DOWNLOAD.out.versions)
        ch_fasta = GENCODE_DOWNLOAD.out.fasta.map { that -> [[id:that.Name], that] }
        ch_gtf = GENCODE_DOWNLOAD.out.gtf.map { that -> [[id:that.Name], that] }
    } else {
        ch_fasta = Channel.fromPath(params.fasta).map { that -> [[id:that.Name], that] }
        ch_gtf = Channel.fromPath(params.gtf).map { that -> [[id:that.Name], that] }
    }

    if (!exists_not_empty(params.fai)){
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = Channel.fromPath(params.fai).map { that -> [[id:that.Name], that] }
    }

    if ((!exists_not_empty(params.hgnc_ref) || !exists_not_empty(params.hgnc_date)) && !params.skip_vcf){
        HGNC_DOWNLOAD( )
        ch_versions = ch_versions.mix(HGNC_DOWNLOAD.out.versions)
        ch_hgnc_ref = HGNC_DOWNLOAD.out.hgnc_ref
        ch_hgnc_date = HGNC_DOWNLOAD.out.hgnc_date
    } else {
        ch_hgnc_ref = Channel.fromPath(params.hgnc_ref).map { that -> [[id:that.Name], that] }
        ch_hgnc_date = Channel.fromPath(params.hgnc_date).map { that -> [[id:that.Name], that] }
    }

    if (!exists_not_empty(params.rrna_intervals)){
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        GET_RRNA_TRANSCRIPTS(ch_gtf)
        ch_versions = ch_versions.mix(GET_RRNA_TRANSCRIPTS.out.versions)
        GATK4_BEDTOINTERVALLIST(GET_RRNA_TRANSCRIPTS.out.bed, GATK4_CREATESEQUENCEDICTIONARY.out.dict )
        ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)
        ch_rrna_interval = GATK4_BEDTOINTERVALLIST.out.interval_list
    } else {
        ch_rrna_interval = Channel.fromPath(params.rrna_intervals).map { that -> [[id:that.Name], that] }
    }

    if (!exists_not_empty(params.refflat)){
        GTF_TO_REFFLAT(ch_gtf)
        ch_versions = ch_versions.mix(GTF_TO_REFFLAT.out.versions)
        ch_refflat = GTF_TO_REFFLAT.out.refflat.map { that -> [[id:that.Name], that] }
    } else {
        ch_refflat = Channel.fromPath(params.refflat).map { that -> [[id:that.Name], that] }
    }

    if (!exists_not_empty(params.salmon_index) || !exists_not_empty(params.salmon_index_stub_check)){ // add condition for qc
        GFFREAD(ch_gtf, ch_fasta.map{ it -> it[1] })
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        SALMON_INDEX(ch_fasta.map{ it -> it[1] }, GFFREAD.out.gffread_fasta.map{ it -> it[1] })
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
        ch_salmon_index = SALMON_INDEX.out.index
    } else {
        ch_salmon_index = Channel.fromPath({params.salmon_index})
    }

    def star_index_tools = tools.intersect(["starindex", "starfusion", "arriba", "ctatsplicing"])
    if (star_index_tools.size() > 0 && (!exists_not_empty(params.starindex_ref) || !exists_not_empty(params.starindex_ref_stub_check))) {
        STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        ch_starindex_ref = STAR_GENOMEGENERATE.out.index
    } else {
        ch_starindex_ref = Channel.fromPath(params.starindex_ref).map { that -> [[id:that.Name], that] }
    }

    if (tools.contains("arriba") && (
            !exists_not_empty(params.arriba_ref_blacklist) ||
            !exists_not_empty(params.arriba_ref_known_fusions) ||
            !exists_not_empty(params.arriba_ref_protein_domains)
        )
    ) {
        ARRIBA_DOWNLOAD(params.genome)
        ch_versions = ch_versions.mix(ARRIBA_DOWNLOAD.out.versions)
        ch_arriba_ref_blacklist = ARRIBA_DOWNLOAD.out.blacklist
        ch_arriba_ref_cytobands = ARRIBA_DOWNLOAD.out.cytobands
        ch_arriba_ref_known_fusions = ARRIBA_DOWNLOAD.out.known_fusions
        ch_arriba_ref_protein_domains = ARRIBA_DOWNLOAD.out.protein_domains
    } else {
        ch_arriba_ref_blacklist = Channel.fromPath(params.arriba_ref_blacklist)
        ch_arriba_ref_cytobands = Channel.fromPath(params.arriba_ref_cytobands)
        ch_arriba_ref_known_fusions = Channel.fromPath(params.arriba_ref_known_fusions)
        ch_arriba_ref_protein_domains = Channel.fromPath(params.arriba_ref_protein_domains)
    }


    if (tools.contains("fusioncatcher") && (!exists_not_empty(params.fusioncatcher_ref) || !exists_not_empty(params.fusioncatcher_ref_stub_check))) {
            FUSIONCATCHER_BUILD(params.genome_gencode_version)
            ch_versions = ch_versions.mix(FUSIONCATCHER_BUILD.out.versions)
            ch_fusioncatcher_ref = FUSIONCATCHER_BUILD.out.reference
    }
    else {
        ch_fusioncatcher_ref = Channel.fromPath(params.fusioncatcher_ref)
    }

    def starfusion_tools = tools.intersect(["starfusion", "ctatsplicing"])
    if (starfusion_tools.size() > 0 && (!exists_not_empty(params.starfusion_ref) || !exists_not_empty(params.starfusion_ref_stub_check))) {
            STARFUSION_BUILD(ch_fasta, ch_gtf, params.fusion_annot_lib, params.species)
            ch_versions = ch_versions.mix(STARFUSION_BUILD.out.versions)
            if (tools.contains("ctatsplicing")) {
                CTATSPLICING_PREPGENOMELIB(
                    STARFUSION_BUILD.out.reference,
                    params.ctatsplicing_cancer_introns
                )
                ch_versions = ch_versions.mix(CTATSPLICING_PREPGENOMELIB.out.versions)
                ch_starfusion_ref = CTATSPLICING_PREPGENOMELIB.out.reference
            } else {
                ch_starfusion_ref = STARFUSION_BUILD.out.reference
            }
    }
    else {
        ch_starfusion_ref = Channel.fromPath(params.starfusion_ref)
    }


    if (tools.contains("fusionreport") && (!exists_not_empty(params.fusionreport_ref) || !exists_not_empty(params.fusionreport_ref_stub_check))) {
        if (!params.no_cosmic && (!params.cosmic_username || !params.cosmic_passwd)) {
            error('COSMIC username and/or password missing, this is needed to download the fusionreport reference')
        }
        FUSIONREPORT_DOWNLOAD()
        ch_versions = ch_versions.mix(FUSIONREPORT_DOWNLOAD.out.versions)
        ch_fusionreport_ref = FUSIONREPORT_DOWNLOAD.out.fusionreport_ref
    } else {
        ch_fusionreport_ref = Channel.fromPath(params.fusionreport_ref).map { that -> [[id:that.Name], that] }
    }

    emit:
    ch_fasta
    ch_gtf
    ch_fai
    ch_hgnc_ref
    ch_hgnc_date
    ch_rrna_interval
    ch_refflat
    ch_salmon_index
    ch_starindex_ref
    ch_arriba_ref_blacklist
    ch_arriba_ref_cytobands
    ch_arriba_ref_known_fusions
    ch_arriba_ref_protein_domains
    ch_fusioncatcher_ref
    ch_starfusion_ref
    ch_fusionreport_ref
    versions        = ch_versions
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def exists_not_empty(path) {
    def path_to_check = file(path)
    return path_to_check.exists() && !path_to_check.isEmpty()
}
