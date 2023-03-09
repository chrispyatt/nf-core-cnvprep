/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowCnvprep.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.bed, params.refGenome, params.map, params.segdup ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.bed) { capture_bed = file(params.bed) } else { exit 1, 'Input capture bed not specified!' }
if (params.refGenome) { ref_genome = file(params.refGenome) } else { exit 1, 'Input reference genome not specified!' }

// Make optional parameters variables
if (params.map) { map_bed = file(params.map) }
if (params.segdup) { segdup_bed = file(params.segdup) }

// Make Groovy map for tuples (may need to change later)
meta_inp = [ id:'test', single_end:false ]

print "\nINPUTS = $ref_genome, $capture_bed, $map_bed\n"



//print ( ref_genome, capture_bed, map_bed )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { GATK4_INDEXFEATUREFILE            } from '../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_PREPROCESSINTERVALS         } from '../modules/nf-core/gatk4/preprocessintervals/main'
include { GATK4_ANNOTATEINTERVALS           } from '../modules/nf-core/gatk4/annotateintervals/main'
include { UNTAR                             } from '../modules/nf-core/untar/main'                                              



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow CNVPREP {

    ch_versions = Channel.empty()


    //
    // SUBWORKFLOW: any local workflow code
    


    //
    // MODULE: Run Untar (on reference genome archive)
    //

    //ref_archive = Channel.of(UNTAR ( [ meta_inp, ref_genome ] ).untar)

    /*
    Channel
    .of(UNTAR ( [ meta_inp, ref_genome ] ).untar)
    .branch { run ->
        fasta: run.toString().endsWith('.fa')
        dict: run.toString().endsWith('.dict')
        fai: run.toString().endsWith('.fai')
    }
    .set { ref_archive }
    */

    untar_out_ch = UNTAR([ meta_inp, ref_genome ] ).untar.map {it ->  it[1]  }
    
    branched_ch = untar_out_ch.branch { run ->
        fasta: run.toString().endsWith('.fa')
        dict: run.toString().endsWith('.dict')
        fai: run.toString().endsWith('.fai')
    }



    print "\nTHIS IS REF_ARCHIVE:\n"
    untar_out_ch.view() { "channel: $it \n" }
    branched_ch.fasta.view() { "fasta: $run \n" }
    //print untar_out_ch
    //print "\nBRANCH FASTA:\n"
    //untar_out_ch.view() { "fasta: $it \n" }


    //print ref_archive.fasta.view()
    //print "\n"
    //print ref_archive.dict.view()
    //print "\n"
    //print ref_archive.fai.view()
    //print "\n"
    //print ref_archive.view()
    //print "\n"

    //
    // MODULE: Run PreprocessIntervals
    //
    
    //prepro_ints = Channel.of(

    /*
    Channel
    .of( GATK4_PREPROCESSINTERVALS ( [ meta_inp, capture_bed ], ref_archive.fasta, ref_archive.dict, ref_archive.fai ) )
    .branch {
        interval_list: it.toString().endsWith('.interval_list')
    }
    .set { prepro_ints }
    //)
    */

    //
    // MODULE: Run IndexFeatureFile
    //

    //to_be_indexed = Channel.of( map_bed, segdup_bed )

    //Channel
    //.of( GATK4_INDEXFEATUREFILE ( [ meta_inp, map_bed ] ) )
    //.set { indexes }


    //.branch {
    //    map_idx: it.toString().endsWith('.idx')[0]
    //    segdup_idx: it.toString().endsWith('.idx')[1]
    //}

    //
    // MODULE: Run AnnotateIntervals
    //
    /*
    anno_ints = GATK4_ANNOTATEINTERVALS (
        [ meta_inp, prepro_ints.interval_list ],
        fasta=ref_archive.fasta,
        dict=ref_archive.dict,
        fai=ref_archive.fai,
        mappable_regions=map_bed,
        mappable_regions_tbi=indexes,
        segmental_duplication_regions='',
        segmental_duplication_regions_tbi=''
        )
    */


}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
