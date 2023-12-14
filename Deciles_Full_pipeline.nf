#!/usr/bin/env nextflow

//*---------------------------------------------
//* Step 1 pre-processing the data 
//*---------------------------------------------
// Calls python script
//Data should be in Peak x Cell type format (peak = rows, cell type = columns) 
//Value can be z-scores/ other metrics such as seurat metadata metrics
//You can set the threshold you want to select for here, it is =/> to the threshold  

all_peak_excel_file = file('/path/to/peakxcell/excel/file')
bedfile_dir = '/path/for/storing/bedfiles'
op0 = '/path/to/Data_prep_threshold_deciles.py'
threshold = '0' 
//pick your own threshold, set to 0 here by default 

process dataPreprocess {
    echo true  
conda '/path/to/conda/env/with/python3.7'
label 'small'
publishDir "$bedfile_dir", pattern: '*.bed', mode: 'copy'


	//Define input 
    input:
    path all_peak_excel_file
    val threshold 


	//Define output 
    output: 
    path "*.bed" into bedfile
    //file "*.bed" into bedfile

script: 
"""
python $op0 $all_peak_excel_file $bedfile_dir $threshold 

"""
}


//*---------------------------------------------------
//* Step 2 liftOver UCSC toold for mapping genomic coordinates
//*---------------------------------------------------

//If you already have .bed files without the need for the first step, you can take the first process out and this pipeline can still be used 
//Just give the input bed files using the below line
//bedfile = Channel.fromPath('/path/to/bedfiles/*.bed')


LP = '/path/for/storing/the/lifted_bed'
bimfile = Channel.fromPath('/path/to/hg38_plink_files/*.bim')
bimpath = '/path/to/hg38_plink_files'
op1 = '/ldsc/make_annot.py'
op2 = '/ldsc/ldsc.py'
// op1 and op2 don't need to be changed if you are using the LDSC singularity image
//if you downloaded LDSC from elsewhere, this needs to be paths to your LDSC
snppath = '/path/to/GRCh38'
baseDir = '/path/to/your/base'
Annot_LD = '/path/to/put/your/annot/and/LD/files'

process liftOver {
    echo true
label 'small'
publishDir "$LP", pattern: 'lifted_*.bed', mode: 'copy'

    input:
		path chain from '/path/to/your/liftOver/chain/file/ex/mm10ToHg38.over.chain'
		//file bedfil 
        //(the above could be used if if .bed files were directly given without step 1)
        file bedfileList from bedfile.flatten()
		val liftover from '/rds/general/user/xyz16/home/liftOver'

    output:
    	path "lifted_*.bed" into liftedbed
  


    script:
    x = bedfileList.name
    """
	$liftover $bedfileList $chain "lifted_${x}" "unlifted_${x}"
    """
    
}

/*---------------------------------------------------
//* Step 3 Create .annot files for S-LDSC
//*--------------------------------------------------
process makeannot {
    echo true
label 'singularity'
label 'small'
publishDir "$Annot_LD", pattern: '*.annot.gz', mode: 'copy'
    input:
        file liftedbed
        each bim from bimfile
        val op1
    output:
        path "*.annot.gz" into annotfile

    script:
    bed = liftedbed.name.replace(".bed", "").replace("lifted_", "")
    chr = bim.name.replace(".bim", "").replace("1000G.EUR.hg38.", "")
//check the names of the chr, because this could be outdaated if plink files update in the future 

    """
    python $op1 --bed-file $liftedbed --bimfile $bim --annot-file "${bed}.${chr}.annot.gz"
    """
}

/*---------------------------------------------------
//* Step 4 Create LD files for S-LDSC
//*--------------------------------------------------

process ldscore {
    echo true
label 'singularity'
//only need the above label if using singularity image 
label 'big'

    publishDir "$Annot_LD", pattern: '*.ldscore.gz', mode: 'copy'
    publishDir "$Annot_LD", pattern: '*.l2.M', mode: 'copy'
    publishDir "$Annot_LD", pattern: '*.l2.M_5_50', mode: 'copy'

    input:
        file annotfile
        val bimpath
        val snppath

    output:
        path "*.ldscore.gz" 
        path "*.l2.M"
        path "*.l2.M_5_50"

    script:
    String[] temp = annotfile.name.split("\\.");
    String bed = temp[0];
    String chr = temp[1];

    """
    python $op2 --l2 --bfile "${bimpath}/1000G.EUR.hg38.${chr}" --ld-wind-cm 1 --annot $annotfile --thin-annot --out "${bed}.${chr}" --print-snps "${snppath}/list.txt"
    """

}

//***-----------------------------------------------------
//*** Step 5 Patrition heritability into different cell types  
//***-----------------------------------------------------

gwas_path = '/path/to/GWAS_sumstats'
gwas_name = 'GWAS.sumstats.gz'

//weight, reference and frq file names 
//should be double checked incase using different versions 

weights_path = '/path/to/hg38_weights'
weight_name = 'weights.hm3_noMHC.'

ref_ld_path = '/path/to/hg38_baselineLD_v2.2'
baseline_name = 'baselineLD.'

frq_path= '/path/to/hg38_plink_files'
frq_name = '1000G.EUR.hg38.'

LDSC_results = '/path/to/store/your/S-LDSC/results'

process patrition_heritability{
    echo true  
    label 'big_mem'
    label 'singularity'
    publishDir "$LDSC_results", pattern: '*.results', mode: 'copy'

    input:
    file liftedbed
    val Annot_LD
    val frq_path
    val frq_name  
    val weights_path 
    val weight_name 
    val gwas_path 
    val gwas_name 
    val ref_ld_path 
    val baseline_name 

    output:
    path "*.results" into results_ch

    script:
    cell_type_name = liftedbed.name.replace(".bed", "").replace("lifted_", "")

    """
    python $op2 --h2 "${gwas_path}/${gwas_name}" --w-ld-chr ${weights_path}/${weight_name} --ref-ld-chr ${Annot_LD}/${cell_type_name}.,${ref_ld_path}/${baseline_name} --overlap-annot --frqfile-chr ${frq_path}/${frq_name} --out ${cell_type_name} --print-coefficient
    """
}

//*---------------------------------------------
//* Step 6 Python function for organizing the LDSC results
//*---------------------------------------------

op3 = '/path/to/home/decile_results_organize.py'
results_table = '/path/to/where/you/want/your/results'

process resultsPreprocess {
    echo true  
conda '/rds/general/user/xyz16/home/anaconda3/envs/py37'
label 'small'
publishDir results_table, mode: 'copy', pattern: 'LDSC_results_deciles.csv'

    //Define input 
    input:
    path LDSC_results


    //Define output 
    output: 
    file("LDSC_results_deciles.csv") into results_table

script: 
"""
python $op3 $LDSC_results

"""
}

//*** --------------------
//* Step 7 R multiple corrections  
//*-----------------------

multiple_corrections_script = '/path/to/NF_multiple_corrections.R'
results_table_final = '/path/to/where/you/want/your/LDSC_summary_table'

process multi_corrections {
conda '/path/to/your/conda/Renv'
// library (FSA) is needed, install if you don't have it 

publishDir results_table_final, mode: 'copy', pattern: 'LDSC_results.csv'
label 'small'

    input:
    file ldscResults from results_table

    output: 
    file("LDSC_results.csv") into results_table_final

    script:
    """
    Rscript $multiple_corrections_script $ldscResults $ldscResults
    """
}

//*** --------------------
//* Step 8 R Plot for S-LDSC results  
//*-----------------------

plotting_script = 'path/to/Plot_LDSC_results_deciles.R'
plot_name = 'LDSC_results_plot_Deciles.pdf'

process plotS-LDSC{
    conda 'path/to/your/conda/Renv'
//libraries include (ggplot2,ggpubr,grid,patchwork,cowplot), ensure they are installed 

publishDir LDSC_plot, mode: 'copy', pattern: 'LDSC_results_plot_Deciles.pdf'
//if you change the plot_name, need to change it in the above line also, and in output

label 'small'

    input: 
    file tableforplot from results_table_final
    val 

    output: 
    file ("LDSC_results_plot_Deciles.pdf") into LDSC_plot 

    script: 
    """
    Rscript $plotting_script $tableforplot $plot_name
    """
}


//*** ----------------------------------- 
//*** Step 9 Examine cCRE overlaps with the ATAC-seq peaks 
//***------------------------------------

cCRE_finder_script = 'path/to/Screen_finder_multiple.py'
cCRE_dir = 'path/to/cCREs'
cCRE_out_path 'path/to/store/overlaps'
//files in the cCRE_dir are best in .txt formats, those used here are downloaded from ENCODE SCREEN project

procress find_cCRE_overlaps {
    conda '/rds/general/user/xyz16/home/anaconda3/envs/py37'
    label 'small'
publishDir "$cCRE_out_path", pattern: '*.csv', mode: 'copy'

 
    input: 
    val cCRE_dir

    output:
    path "*.csv" 

    script: 
    """
    python $cCRE_finder_script $cCRE_dir $cCRE_out_path
    """

// This step could also be taken out and ran independently if you want to test cCRE overlaps with any .bed files
// Bed files should have quantiles attached to them though 
}
















