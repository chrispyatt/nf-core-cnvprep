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

process GET_DX_INPUTS {
    """
    echo "I AM GETTING INPUTS"
    echo ${ref_genome/"dx://"//}
    dx download $ref_genome, $capture_bed, $map_bed
    """
}


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



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow CNVPREP {

    ch_versions = Channel.empty()

    GET_DX_INPUTS()

    //
    // SUBWORKFLOW: any local workflow code
    
    /* this should probably be somewhere else (validate inputs section?)
    process UNPACK {
        """
        tar -xzvf $ref_genome
        """
    }
    UNPACK()
    */

    process TEST {
        """
        ls .
        """
    }
    TEST()

    //
    // MODULE: Run PreprocessIntervals
    //
    /*
    prepro_ints = GATK4_PREPROCESSINTERVALS (
        [ meta_inp, capture_bed ],
        fasta='genome.fa',
        dict='genome.dict',
        fai='genome.fa.fai'
        )
    */

    //
    // MODULE: Run IndexFeatureFile
    //
    /*
    GATK4_INDEXFEATUREFILE (
        [ meta_inp, feature_file=map_bed ]
    )
    */
    /*
    GATK4_INDEXFEATUREFILE (
        feature_file=params.segdup
    )
    */

    //
    // MODULE: Run AnnotateIntervals
    //
    /*
    anno_ints = GATK4_ANNOTATEINTERVALS (
        [ meta_inp, prepro_ints.interval_list ],
        fasta='genome.fa',
        dict='genome.dict',
        fai='genome.fa.fai',
        mappable_regions=map_bed,
        mappable_regions_tbi='map_bed.tbi',
        segmental_duplication_regions=map_bed,
        segmental_duplication_regions_tbi='map_bed.tbi' 
        )
    */

    /*
    process SUCCESS {
        """
        echo "GOT TO THE END FAM. GREAT SUCCESS."
        """
    }
    SUCCESS()
    */

    //
    // MODULE: MultiQC
    //
    /*
    workflow_summary    = WorkflowCnvprep.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCnvprep.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
