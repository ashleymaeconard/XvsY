path="/data/compbio/inathoo/snakemake_test/test7/"

##rule run_deseq2:
##	input: expand(path+"readcounts/metadata_allSamples_full.csv"), expand(path+"readcounts/FB_gene_ByCondition_countTable_full.csv"), expand(path)
##	output: directory(expand(path+"data"))
##	shell:
##		'''
##		Rscript scripts/deseq2_comparison_final.r {input} cRNAi_e 0 0 0.05 dme elavGRFP_e wald 3
##		Rscript scripts/deseq2_comparison_final.r {input} mRNAi_e 0 0 0.05 dme elavGRFP_e wald 3
##		'''

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
		expand(path+"deg_sets/barPlot/sets")
	output:
		directory(expand(path+"deg_sets/gene_beds"))
	shell:
		"./scripts/all_id_to_gene.sh {input}"

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
		"python scripts/global_XAv5.py -f {input} -s {output} -g {path}'scripts'"

rule go_analysis:
	input:
		expand(path+"deg_sets/gene_beds")
	output:
		directory(expand(path+"go_analysis/clusters"))
	shell:
		"./scripts/runGOall.sh {input} {output}"	


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