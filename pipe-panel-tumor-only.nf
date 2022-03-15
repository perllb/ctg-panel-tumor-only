#!/usr/bin/env nextFlow

// PANEL TUMOR ONLY DRAGEN PIPELINE 

// set variables
runfolder = params.runfolder
basedir = params.basedir
metaID = params.metaid
OUTDIR = params.outdir
FQDIR = params.fqdir
QCDIR = params.qcdir
FQCDIR = params.fqcdir
CTGQC = params.ctgqc
demux = params.demux
ffpe = params.ffpe
b2farg = params.bcl2fastq_arg
drarg = params.dragen_arg
index = params.index
nirvanadir = params.nirvanadir

// Read and process sample sheet
sheet = file(params.sheet)

// samplesheet to be parsed as input channel (take everything below [Data])
channel_sheet = file("$basedir/samplesheet.channel.nf.panelTonly.csv")

// create new samplesheet parsed to fit the format for dragen demux
newsheet = "${basedir}/samplesheet.demux.nf.panelTonly.csv"

// Read and process sample sheet
all_lines = sheet.readLines()
write_b = false // if next lines has sample info
channel_sheet.text=""     

for ( line in all_lines ) {

    if ( write_b ) {
	channel_sheet.append(line + "\n")
    }
    if (line.contains("[Data]")) {
	write_b = true
    }
}


println "========================================="
println ">>> panel T only somatic dragen pipeline "
println ""
println "> INPUT: "
println "> runfolder		: $runfolder "
println "> sample-sheet		: $sheet "
println "> run-meta-id		: $metaID "
println "> basedir		: $basedir "
println "> bcl2fastq args	: $b2farg "
println ""
println "> OUTPUT: "
println "> output-dir		: $OUTDIR "
println "> fastq-dir		: $FQDIR "
println "> qc-dir		: $QCDIR "
println "> fastqc-dir		: $FQCDIR "
println "> ctg-qc-dir		: $CTGQC "
println "============================="

// sample info
Channel
    .fromPath(channel_sheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Ref, row.panel, row.annotate, row.fastq_1, row.fastq_2, row.bam) }
    .unique()
    .tap{infoSamples}
    .into{ move_fastq_csv; analyze_csv; fastqc_go }

// project info
Channel
    .fromPath(channel_sheet)
    .splitCsv(header:true)
    .map { row -> tuple(  row.Sample_Project, row.panel ) }
    .unique()
    .tap{infoProjects}
    .set{ dragen_summary  }

println " > Samples to process: "
infoSamples.subscribe{ println "Sample: $it" }


println " > Projects to process: "
infoProjects.subscribe{ println "Sample: $it" }

// Parse samplesheet
process parsesheet {

	tag "$metaID"

	input:
	val "start"

	output:
	val newsheet into demux_sheet

	when:
	demux == 'y'

	"""
python $basedir/bin/ctg-parse-samplesheet.dragen-panel.py -s $sheet -o $newsheet -i $index
	"""
}

// dragen demux
process demux {

    tag "$metaID"	
    label 'dragen'

    input:
    val newsheet from demux_sheet

    output:
    val "x" into mv_fastq

    when:
    demux = "y"
        
    """
    export LC_ALL=C        

    mkdir -p ${FQDIR}
    /opt/edico/bin/dragen --force --bcl-conversion-only=true \\
           --bcl-input-directory ${runfolder} \\
	   --output-directory ${FQDIR} \\
	   --sample-sheet ${newsheet} \\
	   --no-lane-splitting true \\
	   ${b2farg}

     """
}

process moveFastq {

    tag "${projid}_${sid}"

    input:
    val x from mv_fastq
    set sid, projid, ref, panel, annotate, fastq_1, fastq_2, bam from move_fastq_csv

    output:
    val "y" into run_analysis
    val "x" into run_qc

    when:
    demux = 'y'

    """
    mkdir -p ${OUTDIR}/${projid}
    mkdir -p ${OUTDIR}/${projid}/fastq
   
    # If there is a directory per project
    if [ -d ${FQDIR}/${projid}/ ]; then
       mv ${FQDIR}/${projid}/${fastq_1} ${OUTDIR}/${projid}/fastq/
       mv ${FQDIR}/${projid}/${fastq_2} ${OUTDIR}/${projid}/fastq/
    else
       mv ${FQDIR}/${fastq_1} ${OUTDIR}/${projid}/fastq/
       mv ${FQDIR}/${fastq_2} ${OUTDIR}/${projid}/fastq/
    fi     

    """

}

// Channel to start analysis if demux == 'n'
// Projects
if ( demux == 'n' ) {
   Channel
	 .value("1")
    	 .into{ run_analysis; run_qc }
}

// dragen run : align, vc + metrics
process dragen_align_vc {

    tag "${projid}_${sid}"
    label 'dragen' 

    input:
    val x from run_analysis
    set sid, projid, ref, panel, annotate, fastq_1, fastq_2, bam from analyze_csv

    output:
    val x into done_analyze
    val projid into dragen_metrics
    set sid, projid, ref, annotate, val("${OUTDIR}/${projid}/dragen/vcf/${sid}/${sid}.hard-filtered.vcf.gz"), val("${OUTDIR}/${projid}/dragen/vcf/${sid}/tumorSV.vcf.gz") into annotate_vcf

    """
    export LC_ALL=C

    # Get target panel file
    if [ $panel == "comprehensive" ] || [ $panel == "Comprehensive" ]
    then
	if [ $ref == 'hg38' ]; then
	    targetfile=${params.target_twist_comprehensive_hg38}
	elif [ $ref == 'mm10' ]; then
	    targetfile=${params.target_twist_comprehensive_mm10}
        fi
    elif [ $panel == "core" ] || [ $panel == "Core" ]
    then
        if [ $ref == 'hg38' ]; then
	    targetfile=${params.target_twist_core_hg38}
	elif [ $ref == 'mm10' ]; then
	    targetfile=${params.target_twist_core_mm10}
        fi
    elif [ $panel == "gmck" ] || [ $panel == "GMCK" ]
    then
        if [ $ref == 'hg19' ]; then
	    targetfile=${params.target_gmck_hg19}
	else
	    echo 'ERROR: GMCK panel can only be run with hg19 reference genome'
        fi
    elif [ $panel == "gms570" ] || [ $panel == "GMS570" ]
    then
        if [ $ref == 'hg19' ]; then
	    targetfile=${params.target_gms570_hg19}
	else
	    echo 'ERROR: GMS570 panel can only be run with hg19 reference genome'
        fi
    elif [ $panel == "custom"  ] || [ $panel == "Custom" ] 
    then
        targetfile=${params.custom_target}
    else
        echo '>PANEL NOT RECOGNIZED!'
	echo 'in samplesheet - only 'comprehensive', 'core' and 'custom' can be specified in 'panel' section'
        targetfile='ERR'
    fi

    # Get species based on ref
    if [[ $ref == hg* ]]
    then 
    	species='human'
    elif [[ $ref == mm* ]]
    then
	species='mouse'
    else
        echo 'species cannot be extracted from reference: $ref'
	exit 1;
    fi

    dragendir=${OUTDIR}/${projid}/dragen/
    mkdir -p \$dragendir
    mkdir -p \${dragendir}/metrics
    mkdir -p \${dragendir}/vcf/
    mkdir -p \${dragendir}/vcf/${sid}

    outdir=\${dragendir}/metrics/${sid}
    mkdir -p \$outdir

    R1=${OUTDIR}/${projid}/fastq/$fastq_1
    R2=${OUTDIR}/${projid}/fastq/$fastq_2

    echo "R1: '\${R1}'"
    echo "R2: '\${R2}'"
    echo "sid: '${sid}'"
    echo "padding: '${params.padding}'"
    echo "outdir: '\${outdir}'"
    echo "targetfile: '\${targetfile}'"
    echo "species: '\${species}'"   
    echo "ffpe: '´{params.ffpe}'"
   
    /opt/edico/bin/dragen -f -r /staging/\$species/reference/$ref \\
        --tumor-fastq1 \${R1} \\
        --tumor-fastq2 \${R2} \\
        --RGID-tumor ${projid}_${sid} \\
        --RGSM-tumor $sid \\
        --intermediate-results-dir /staging/tmp/ \\
        --enable-map-align true \\
        --enable-map-align-output true \\
        --vc-target-bed \$targetfile \\
        --vc-target-bed-padding ${params.padding} \\
        --output-format bam \\
        --output-directory \$outdir \\
        --enable-variant-caller true \\
        --enable-sv true \\
        --output-file-prefix $sid \\
        --remove-duplicates true \\
        --qc-coverage-region-1 \$targetfile \\
	--qc-coverage-region-padding-1 ${params.padding} \\
        --qc-coverage-ignore-overlaps true \\
        --qc-coverage-filters-1 "mapq<20,bq<20" \\
	--vc-enable-orientation-bias-filter ${params.ffpe} \\
	--vc-sq-filter-threshold 20 --vc-enable-af-filter true \\
	${drarg}
	
    # move vcfs to vcf dir     
    # SNV
    mv \${outdir}/${sid}*.vcf.gz* \${dragendir}/vcf/${sid}/
    # SV
    mv \${outdir}/sv/results/variants/*.vcf.gz* \${dragendir}/vcf/${sid}/

    #move bams to bams dir
    mkdir -p \${dragendir}/bams
    mkdir -p \${dragendir}/bams/${sid}
    mv \${outdir}/${sid}*.bam* \${dragendir}/bams/${sid}/
    """
}

// Annotate vcf
process annotate {

	tag "${projid}-${sid}"
	label 'dragen' 	

	input:
	set sid, projid, ref, annotate, snv, sv from annotate_vcf
	
	output:
	set sid, projid, ref, annotate into filter_vcf

	when:
	annotate == 'y'

	"""
	nirvanaref=$nirvanadir/$ref/

	vcfdir='${OUTDIR}/${projid}/dragen/vcf/$sid'

	outSNV_BN=\$(basename $snv .vcf.gz)
	outSV_BN=\$(basename $sv .vcf.gz)
	
	outSNV=\$vcfdir/\${outSNV_BN}.annotated.nirvana.txt
	outSV=\$vcfdir/\${outSV_BN}.annotated.nirvana.txt

	# Convert ref 
	# hg38 to GRCh38 
	if [ $ref == 'hg38' ]; then
		Gref='GRCh38'
	fi		
	
	# hg19 to GRCh37 
	if [ $ref == 'hg19' ]; then
		Gref='GRCh37'
	fi		

	# SNV annotate
	/opt/edico/share/nirvana/Nirvana -c \$nirvanaref/Cache/\$Gref/Both -r \$nirvanaref/References/Homo_sapiens.\$Gref.Nirvana.dat --sd \$nirvanaref/SupplementaryAnnotation/\$Gref -i $snv -o \$outSNV

	# SV annotate		 
	/opt/edico/share/nirvana/Nirvana -c \$nirvanaref/Cache/\$Gref/Both -r \$nirvanaref/References/Homo_sapiens.\$Gref.Nirvana.dat --sd \$nirvanaref/SupplementaryAnnotation/\$Gref -i $sv -o \$outSV

	"""
}


// fastqc 
process fastqc {

	tag "${projid}_${sid}"
	
	input:
	val x from run_qc
	set sid, projid, ref, panel, annotate, fastq_1, fastq_2, bam from fastqc_go

        output:
        val projid into multiqc_fastqc

	
	"""

        qcdir=${OUTDIR}/${projid}/qc
        fqcdir=${OUTDIR}/${projid}/qc/fastqc
        mkdir -p \$qcdir
        mkdir -p \$fqcdir 

	R1=${OUTDIR}/${projid}/fastq/$fastq_1
	R2=${OUTDIR}/${projid}/fastq/$fastq_2

	fastqc -t $task.cpus --outdir \$fqcdir \$R1
	fastqc -t $task.cpus --outdir \$fqcdir \$R2
 	
	"""
    
}

process dragen_stats {

        tag "${projid}"

	input: 
	val "y" from dragen_metrics.collect()
	set projid, panel from dragen_summary.unique()
	
	output:
	val projid into multiqc_dragen

	"""

    if [ $panel == "comprehensive" ] || [ $panel == "Comprehensive" ]
    then
        targetfile='TWIST-comprehensive'
    elif [ $panel == "core" ] || [ $panel == "Core" ]
    then
        targetfile='TWIST-core'
    else
        targetfile='Custom panel'
    fi

	mkdir -p ${OUTDIR}/$projid/qc/dragen
	${basedir}/bin/ctg-dragen-stats-panel -p $projid -i ${OUTDIR}/$projid/dragen/ -o ${OUTDIR}/$projid/qc/dragen/ -a $params.padding -t \${targetfile}
	
	"""
}



process multiqc {

    tag "${projid}"

    input:
    set projid, projid2 from multiqc_fastqc.unique().phase(multiqc_dragen.unique())

    output:
    val projid into pipedone

    """
    
    cd $OUTDIR
    multiqc -f ${OUTDIR}/$projid/ --outdir ${OUTDIR}/$projid/qc/multiqc/ -n ${projid}_panelTonlySomatic_dragen_report.html

    mkdir -p ${CTGQC}/$projid

    cp -r ${OUTDIR}/$projid/qc ${CTGQC}/$projid/

    """

}


// write to cronlog when pipeline is ready
process pipe_done {

        tag "$projid"	

	input:
	val projid from pipedone

	""" 

	touch $runfolder/ctg.panel-tumor-only.$projid.done

	cronlog="/projects/fs1/shared/ctg-cron/ctg-pipe-cron/logs/panel-tumor-only/cron-panel-tumor-only.log"
	cronlog_all="/projects/fs1/shared/ctg-cron/ctg-cron.log"
	
	rf=\$(basename $runfolder)
  	echo "\$(date): \$rf: DONE: panel-tumor-only ($projid)" >> \$cronlog
    	echo "\$(date): \$rf: DONE: panel-tumor-only ($projid)" >> \$cronlog_all

	"""

}
