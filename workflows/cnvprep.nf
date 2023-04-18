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


    //
    // MODULE: Run Untar (on reference genome archive)
    //

    untar_out_ch = UNTAR([ meta_inp, ref_genome ] ).untar.multiMap { it ->  
        fasta: it[1][1]
        dict: it[1][0]
        fai: it[1][2] 
        }

    
    fasta_ch = untar_out_ch.fasta
    dict_ch = untar_out_ch.dict
    fai_ch = untar_out_ch.fai

    fasta_ch.view() { "fasta: $it \n" }
    dict_ch.view() { "dict: $it \n" }
    fai_ch.view() { "fai: $it \n" }
    

    //
    // MODULE: Run PreprocessIntervals
    //
    
    prepro_ints = GATK4_PREPROCESSINTERVALS ( [ meta_inp, capture_bed ], fasta_ch, dict_ch, fai_ch ).interval_list.multiMap { it -> interval: it[1] }

    interval_ch = prepro_ints.interval

    interval_ch.view() { "interval_ch: $it \n" }

    /*
    Channel
    .of( GATK4_PREPROCESSINTERVALS ( [ meta_inp, capture_bed ], untar_out_ch.fasta, untar_out_ch.dict, untar_out_ch.fai ) )
    .branch {
        interval_list: it.toString().endsWith('.interval_list')
    }
    .set { prepro_ints }
    */

    prepro_ints.view() { "prepro: $it \n" }

    untar_out_ch.fasta.view() { "fasta: $it \n" }
    untar_out_ch.dict.view() { "dict: $it \n" }
    untar_out_ch.fai.view() { "fai: $it \n" }

    //
    // MODULE: Run IndexFeatureFile
    //

    to_be_indexed_ch = Channel.of( map_bed )

    indexes = GATK4_INDEXFEATUREFILE ( [ meta_inp, map_bed ] ).index.map { it -> it[1] }

    /*
    Channel
    .of( GATK4_INDEXFEATUREFILE ( [ meta_inp, map_bed ] ) )
    .set { indexes }
    */

    indexes.view() { "indexes: $it \n" }


    //.branch {
    //    map_idx: it.toString().endsWith('.idx')[0]
    //    segdup_idx: it.toString().endsWith('.idx')[1]
    //}

    //
    // MODULE: Run AnnotateIntervals
    //
    /*
    fasta_ch.view() { "fasta: $it \n" }
    dict_ch.view() { "dict: $it \n" }
    fai_ch.view() { "fai: $it \n" }
    
    anno_ints = GATK4_ANNOTATEINTERVALS (
        [ meta_inp, prepro_ints ],
        fasta_ch,
        dict_ch,
        fai_ch,
        map_bed,
        indexes,
        '',
        ''
        )
    */
    //anno_ints.view() //{ "annotations: $it \n" }

    process publishOutputs {
        publishDir '.'
        script:
        '''
        echo "this is some test output" > ./test_output.txt
        '''
    }

    //publishOutputs()


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
