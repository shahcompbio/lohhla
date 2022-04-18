# Shared configuration
configfile: "config/v6.yaml"

# Sub-pipelines
include: "rules/common.smk"
include: "rules/hla_loh.smk"

rule all:
	input:
		expand(
			'{outdir}/{subdir}/hla_loh/{{patient_id}}'.format(
				outdir=config['outdir']['hla-loh'],
				subdir=config['outputs']['out'],
			),
			patient_id=config['patients']
		),
