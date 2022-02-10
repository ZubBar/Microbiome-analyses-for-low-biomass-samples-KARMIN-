# Microbiome-analyses-for-low-biomass-samples-KARMIN-

#   Sequence-based filtering pipeline for low biomass samples

Filter the sequences of your blank control(s) from your demultiplexed fasta file.
 
    filter_fasta.py 
	    -f $PWD/seqs.fna 
	    -p Blank 
	    -o $PWD/Blank_seqs.fna
	
Cluster the sequences from the Blank control(s) at 99% 

    pick_otus.py 
	    -i $PWD/Blank_seqs.fna 
	    -m uclust 
	    -o $PWD/otus_Blank.only_uclust_99 
	    -s 0.99 
	    -g 2 
	    --threads 8
	
Use the most representative sequences from the largest clusters (i.e. the “X” top clusters witch contain >1% of the total reads).  

    awk '$0=$1": "NF' $PWD/otus_Blank.only_uclust_99/Blank_seqs_otus.txt |sort -k 2 -n |tail -n X >$PWD/Blank_Xlargestclusters.txt

    pick_rep_set.py 
	    -i $PWD/otus_Blank.only_uclust_99/Blank_seqs_otus.txt 	
	    -f $PWD/Blank_seqs.fna 
	    -m most_abundant 
	    -o $PWD/rep_set_Blank_only.fna
	
    awk '{print $1}' $PWD/Blank_Xlargestclusters.txt |tr -s ':' ' ' |tr -d "[:blank:]" >$PWD/Blank_Xlargestclusters.txt_IDs.txt

    filter_fasta.py 
	    -f $PWD/rep_set_Blank_only.fna 
	    -o $PWD/rep_set_Blank_Xlargestclusters.fna 
	    -s $PWD/Blank_Xlargestclusters.txt_IDs.txt

    --> At this point you have picked the most representative sequences assigned to the X largest clusters of your blank control(s).	

Obtain all reads from these representative sequences in these X clusters.

    grep -Ff $PWD/Blank_Xlargestclusters.txt_IDs.txt $PWD/otus_Blank.only_uclust_99/Blank_seqs_otus.txt | awk '$1=" "' |awk 'BEGIN {OFS="\n"}{$1=$1;print}' >$PWD/Blank_Xlargestclusters_readsIDs.txt

    filter_fasta.py 
	    -f $PWD/Blank_seqs.fna 
	    -o $PWD/Blank_seqs_XlargestclustersAll.fna 
	    -s $PWD/Blank_Xlargestclusters_readsIDs.txt

    grep '^[^>]' $PWD/Blank_seqs_XlargestclustersAll.fna |sort|uniq >$PWD/Blank_seqs_XlargestclustersAll_uniq.fna

Filter now your samples from your demultiplexed fasta file.

    filter_fasta.py 
	    -f $PWD/seqs.fna 
	    -p Samples 
	    -o $PWD/Samples_seqs.fna


Filter your Blank "selected sequences" from your samples 

    grep -vFf $PWD/Blank_seqs_XlargestclustersAll_uniq.fna $PWD/Samples_seqs.fna >$PWD/Samples.Blanklargestclusters_fil.fna

    grep '^[^>]' -B 1 $PWD/Samples.Blanklargestclusters_fil.fna |sed '/--/d' >$PWD/Samples.BlanklargestCl_filt_seqsFinal.fna

    --> with this last fasta file you can now keep on perfomrmig your taxonomic analyses with your standard pipeline.


# Qiime v1.9.1 scripts workflow from http://qiime.org/scripts/

Analyses of taxonomic distribution

Use the open-reference OTU picking command for sequence clustering and OTU assigment (if you aim to use any other database for taxonomy assigment define a parameter file, see http://www.qiime.org/documentation/file_formats.html#qiime-parameters)

 	pick_open_reference_otus.py 
		-i $PWD/seqs.fna 
		-p params.txt
		-o $PWD/pick_open_ref &
	
Retain those OTUs that show a minimum total observation count of 10 (n=10)

	filter_otus_from_otu_table.py 
		-i $PWD/otu_table_mc2_w_tax_no_pynast_failures.biom 
		-n 10
		-o $PWD/otu_table_n10.biom 

Estimate the fraction of the total observation (sequence) count to apply as the minimum total observation count of an OTU to be retained. Eg.: 0.2% 

	filter_otus_from_otu_table.py 
		-i $PWD/otu_table_n10.biom  
		--min_count_fraction 0.002
		-o $PWD/otu_table_n10_mFr02.biom
	
If you aim to group your samples regarding any variable of your study: collapse your samples
	
	collapse_samples.py 
		-b $PWD/otu_table_n10_mFr02.biom 
		-m $PWD/Mapping_file.txt 
		--output_biom_fp $PWD/collapsed_normed.biom 
		--output_mapping_fp $PWD/collapsed_map.txt .
		--collapse_mode "choose your mode (median/mean)" 
		--collapse_fields "choose your variable" 
		--normalize
	
Bacterial taxonomic distribution bar plots (if you would like to have your samples sorted use option -s)

	summarize_taxa_through_plots.py 
		-i $PWD/collapsed_normed.biom 
		-o $PWD/taxa 
		-s
	
Filter samples from your data set

	filter_samples_from_otu_table.py 
		-i $PWD/otu_table_n10_mFr02.biom 
		-m $PWD/Mapping_file.txt 
		-s 'define your variable'
		-o $PWD/otu_table_n10_mFr02_Variable.biom 
	
alpha_diversity 

	alpha_rarefaction.py 
		-i $PWD/otu_table_n10_mFr02.biom
		-m $PWD/Mapping_file.txt   
		-t rep_set.tre &
		-e "define the upper limit of rarefraction depths"
		-o alpha_diversity
	
beta_diversity

	beta_diversity_through_plots.py 
		-i $PWD/otu_table_n10_mFr02.biom     
		-t rep_set.tre
		-m $PWD/Mapping_file.txt
		-e "define sequencing depth"
		-o $PWD/beta_diversity
	
make 2d plots (optional)

	make_2d_plots.py 
		-i $PWD/beta_diversity/(un)weighted_unifrac_pc.txt 
		-m $PWD/Mapping_file.txt
		-b 'define your variable if desired' 
		-o $PWD/2d_plots 
 
perform a statistical test (optional)

	compare_categories.py 
		-i $PWD/beta_diversity/(un)weighted_unifrac_dm.txt 
		-m $PWD/Mapping_file.txt
		-c 'define your variable' -o permanova_stats 
		--method "choose your statistical test:permanova/anosim/adonis/permdisp/db-RDA...see http://qiime.org/scripts/compare_categories.html)
		-o $PWD/statistic
	
	(suggestion: perform permdisp statistical test before permanova to test for group dispersion homogeneity)

Compare OTU frequencies across sample groups (statistical test eg.: kruskal wallis)

	group_significance.py 
		-i $PWD/otu_table_n10_mFr02_L6.biom 
		-m $PWD/Mapping_file.txt 
		-c "define your variable" 
		-s "choose your statistical test" 
		-o $PWD/statistical_test.txt
  
	For more detailed information about this pipeline referred to "Caporaso, J., Kuczynski, J., Stombaugh, J. et al. QIIME allows analysis of high-throughput community sequencing data. Nat Methods 7, 335–336 (2010). https://doi.org/10.1038/nmeth.f.303"

# SparCC pipeline adaptation from Jonathan Friedman (https://github.com/JCSzamosi/SparCC3) 

Correlation between all OTUs_quantification

	python SparCC.py $PWD/taxa/otu_table_n10_mFr02_L7.txt -c sparcc_cor.txt -v sparcc_cov.txt

Pseudo p-value calculation via bootstrap 

	python MakeBootstraps.py $PWD/taxa/otu_table_n10_mFr02_L7.txt -p simulated_test/ -t permuted_#

	mkdir boot_corr

	mkdir boot_cov

	for i in `seq 0 99`; do python SparCC.py simulated/permuted_$i -c boot_corr/simulated_sparcc_$i.txt -v boot_cov/simulated_sparcc_$i.txt >> boot_sparcc.log; done

	python PseudoPvals.py sparcc_cor.txt boot_corr/simulated_sparcc_#.txt 100 -o $PWD/pvalues/one_sided.txt -t one_sided (or two_sided)

	For more detailed information about these calculations referred to "Friedman, J., & Alm, E. J. (2012). Inferring correlation networks from genomic survey data."

Replace $PWD with your path

Give your option wherever you find ""
