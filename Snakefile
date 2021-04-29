configfile: "config.yaml"

path=config['outdir'] #should end with a '/'

# rule run_deseq2_2ex:
# 	# input: expand(path+"{metadata}", metadata=config["metadata"]), expand(path+"{counts}", counts=config["counts"]), expand(path)
# 	# input: expand(path)
# 	output: directory(expand(path+"data"))
# 	shell:
# 		'''
# 		Rscript scripts/deseq2_comparison_final.r {config[metadata]} {config[counts]} {path} {config[condition]} {config[batch_effect]} {config[time_course]} {config[pval_threshold]} {config[organism]} {config[control]} {config[stat_test]} {config[read_threshold]}
# 		Rscript scripts/deseq2_comparison_final.r {config[metadata]} {config[counts]} {path} {config[condition2]} {config[batch_effect]} {config[time_course]} {config[pval_threshold]} {config[organism]} {config[control2]} {config[stat_test]} {config[read_threshold]}
# 		'''

rule run_deseq2_3ex:
	output: directory(expand(path+"data"))
	shell:
		'''
		Rscript scripts/deseq2_comparison_final.r {config[metadata]} {config[counts]} {path} {config[condition]} {config[batch_effect]} {config[time_course]} {config[pval_threshold]} {config[organism]} {config[control]} {config[stat_test]} {config[read_threshold]}
		Rscript scripts/deseq2_comparison_final.r {config[metadata]} {config[counts]} {path} {config[condition2]} {config[batch_effect]} {config[time_course]} {config[pval_threshold]} {config[organism]} {config[control2]} {config[stat_test]} {config[read_threshold]}
		Rscript scripts/deseq2_comparison_final.r {config[metadata]} {config[counts]} {path} {config[condition3]} {config[batch_effect]} {config[time_course]} {config[pval_threshold]} {config[organism]} {config[control3]} {config[stat_test]} {config[read_threshold]}
		'''

rule get_ids:
	input:
		expand(path+"data")
	output:
		directory(expand(path+"deg_sets/gene_ids"))
	shell:
		"./scripts/all_csv_to_id.sh {input}"

rule find_intersections:
	input:
		expand(path+"deg_sets/gene_ids")
	output:
		directory(expand(path+"deg_sets/barPlot/sets"))
	shell:
		"intervene upset -i {input}/*.txt --type list --save-overlaps --figtype png -o '{input}/../barPlot'"

rule get_gene:
	input:
		# expand(path+"deg_sets/barPlot/sets"), expand(path+"{gtf_path}", gtf_path=config["gtf_path"])
		expand(path+"deg_sets/barPlot/sets")
	output:
		directory(expand(path+"deg_sets/gene_beds"))
	shell:
		"./scripts/all_id_to_gene.sh {input} {config[gtf_path]}"

rule make_boxplots:
	input:
		expand(path+"deg_sets/gene_beds")
	output:
		directory(expand(path+"comparisons/overall_boxplots"))
	shell:
		"python scripts/allBoxPlots.py -f {input} -d {path}'data' -s {output}"	

rule make_XAboxplots:
	input:
		expand(path+"deg_sets/gene_beds")
	output:
		directory(expand(path+"comparisons/xa_boxplots"))
	shell:
		"python scripts/xaBoxPlotsfinal.py -f {input} -d {path}'data' -s {output}"	

rule global_boxplots:
	input:
		expand(path+"data")
	output:
		directory(expand(path+"comparisons/global_boxplots"))
	shell:
		"python scripts/global_XAv5.py -f {input} -s {output} -g {config[gtf_bed]}"

rule go_analysis:
	input:
		expand(path+"deg_sets/gene_beds")
	output:
		directory(expand(path+"go_analysis/clusters"))
	shell:
		"./scripts/runGOall.sh {input} {output} {config[sep_tps]} {config[organism]} {config[go_pval]}"	

rule go_summary:
	input:
		expand(path+"go_analysis/clusters")
	output:
		directory(expand(path+"go_analysis/summary"))
	shell:
		'''
		python scripts/goList.py -f {input} -g BP
		python scripts/goList.py -f {input} -g CC
		python scripts/goList.py -f {input} -g MF
		'''

rule meme_suite_prep:
	input:
		expand(path+"deg_sets/gene_beds")
	output:
		directory(expand(path+"motif_analysis/clusters"))
	shell:
		"./scripts/memesuitePrep.sh {input} {output} {config[reform_genes]} {config[chrom_fa]} {config[tss_only]} {config[organism]}"

rule run_meme:
	input:
		expand(path+"motif_analysis/clusters")
	output:
		directory(expand(path+"motif_analysis/clusters_meme"))
	shell:
		'''
		cp -r {input} {input}_meme
		./scripts/runMeme.sh {input}_meme
		'''

rule run_fimo:
	input:
		expand(path+"motif_analysis/clusters")
	output:
		directory(expand(path+"motif_analysis/clusters_fimo"))
	shell:
		'''
		cp -r {input} {input}_fimo
		./scripts/runFimo.sh {input}_fimo {config[pwm_path]}
		python scripts/fimo_summary.py -p {output}
		'''
