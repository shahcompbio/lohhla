workflow_name = 'hla-loh'

output_dir = config['outdir'][workflow_name]
workspace_dir = config['workspaces'][workflow_name]

def get_wgs_input_paths(wildcards):
    return {
        'normal_bam': wgsruninfo.paths[(wgsruninfo.paired_dict[(wildcards.patient_id, 'NORMAL')], 'WGS-ALIGNMENT', 'bam')],
        'normal_bai': wgsruninfo.paths[(wgsruninfo.paired_dict[(wildcards.patient_id, 'NORMAL')], 'WGS-ALIGNMENT', 'bai')],
        'tumor_bam': wgsruninfo.paths[(wgsruninfo.paired_dict[(wildcards.patient_id, 'TUMOR')], 'WGS-ALIGNMENT', 'bam')],
        'tumor_bai': wgsruninfo.paths[(wgsruninfo.paired_dict[(wildcards.patient_id, 'TUMOR')], 'WGS-ALIGNMENT', 'bai')],
        'alleles': wgsruninfo.paths[(wgsruninfo.paired_dict[(wildcards.patient_id, 'NORMAL')], 'WGS-POLYSOLVER', 'winners_hla')],
    }

rule run_lohhla:
	input:
		unpack(get_wgs_input_paths)
	output:
		# cna="{outdir}/{subdir}/hla_loh/{{patient_id}}/copy_number_loc.txt".format(
		# 	outdir=output_dir,
		# 	subdir=config['outputs']['out']
		# ),
		dir=directory(
			'{outdir}/{subdir}/hla_loh/{{patient_id}}'.format(
				outdir=output_dir,
				subdir=config['outputs']['out']
			),
		),
	params:
		name='lohhla-patient-{patient_id}',
		patient_id='{patient_id}',
		hla_dat=config['hla_loh']['lohhla']['hla_dat'],
		hla_fasta=config['hla_loh']['lohhla']['hla_fasta'],
		gatk_dir=config['hla_loh']['lohhla']['gatk_dir'],
		picard_dir=config['hla_loh']['lohhla']['picard_dir'],
		novo_dir=config['hla_loh']['lohhla']['novo_dir'],
		poly_dir=config['hla_loh']['lohhla']['polysolver_dir'],
		caller=config['hla_loh']['lohhla']['caller'],
		mapping_step=config['hla_loh']['lohhla']['mapping_step'],
		min_coverage_filter=config['hla_loh']['lohhla']['min_coverage_filter'],
		fishing_step=config['hla_loh']['lohhla']['fishing_step'],
		cleanup=config['hla_loh']['lohhla']['cleanup'],
		workspace=workspace_dir,
	log:
		'{outdir}/{subdir}/hla_loh/{{patient_id}}.log'.format(
			outdir=output_dir,
			subdir=config['outputs']['log']
		),
	benchmark:
		'{outdir}/{subdir}/hla_loh/{{patient_id}}.txt'.format(
			outdir=output_dir,
			subdir=config['outputs']['bench']
		),
	shadow:
		'shallow'
	singularity:
		# 'docker://cmopipeline/lohhla:1.1.6'
		'docker://mleventhal/lohhla:latest'
	shell:
		'export PATH=$PATH:/juno/work/shah/vazquezi/software/novocraft; '
		'mkdir -p bam; '
		'mkdir -p flagstat; '
		'mkdir -p {output.dir}; '
		'ln -sf {input.normal_bam} bam/normal.bam &> {log}; '
		'ln -sf {input.normal_bai} bam/normal.bam.bai &>> {log}; '
		'ln -sf {input.tumor_bam} bam/tumor.bam &>> {log}; '
		'ln -sf {input.tumor_bai} bam/tumor.bam.bai &>> {log}; '
		'Rscript {params.workspace}/R/LOHHLAscript.R '
		'--patientId {params.patient_id} '
		'--outputDir {output.dir} '
		'--normalBAMfile bam/normal.bam '
		'--BAMDir bam/ '
		'--hlaPath {input.alleles} '
		'--HLAfastaLoc {params.hla_fasta} '
		'--novoDir {params.novo_dir} '
		'--gatkDir {params.gatk_dir} '
		'--HLAexonLoc {params.hla_dat} '
		'--mappingStep {params.mapping_step} '
		'--minCoverageFilter {params.min_coverage_filter} '
		'--fishingStep {params.fishing_step} '
		'--cleanUp {params.cleanup} '
		# '--CopyNumLoc {output.cna} '
		'&>> {log}; '
