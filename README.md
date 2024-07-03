# **Single-cell deconvolution methods - practical session**  

***

June 28th 2024  

-Datasets and associated metadata: PDAC TCGA (bulk RNA-seq downsampled, PDAC single-cell reference)

-Inference of cell-types using BayesPrism, Bisque, CIBERSORTx, and EPIC 

-Open the sc_deconvolution.Rmd file to continue with the code. You have all the necessary R libraries and input data in your virtual machine.

***

## Resources/Vignettes:  

Bisque vignette: http://127.0.0.1:27074/library/BisqueRNA/doc/bisque.html  
BayesPrism tutorial: https://www.bayesprism.org/pages/tutorial_deconvolution  
MuSic tutorial: https://xuranw.github.io/MuSiC/articles/MuSiC.html  
Omnideconv mouse https://omnideconv.org/immunedeconv/articles/detailed_example_mouse.html  
Omnideconv: https://github.com/omnideconv/omnideconv
CIBERSORTx: https://cibersortx.stanford.edu/  
EcoTyper: https://ecotyper.stanford.edu/carcinoma/
SimBu simulator: https://github.com/omnideconv/SimBu



+ Avila Cobos, F., Alquicira-Hernandez, J., Powell, J.E. et al. Benchmarking of cell type deconvolution pipelines for transcriptomics data. Nat Commun 11, 5650 (2020). https://doi.org/10.1038/s41467-020-19015-1  

+ Huuki-Myers, L. A., Montgomery, K. D., Kwon, S. H., Cinquemani, S., Eagles, N. J., Gonzalez-Padilla, D., Maden, S. K., Kleinman, J. E., Hyde, T. M., Hicks, S. C., Maynard, K. R., & Collado-Torres, L. (2024). Benchmark of cellular deconvolution methods using a multi-assay reference dataset from postmortem human prefrontal cortex. bioRxiv : the preprint server for biology, 2024.02.09.579665. https://doi.org/10.1101/2024.02.09.579665

+ Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., List, M., & Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics (Oxford, England), 35(14), i436–i445. https://doi.org/10.1093/bioinformatics/btz363  

+ Jew, B., Alvarez, M., Rahmani, E. et al. Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6  
  
+ Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022). https://doi.org/10.1038/s43018-022-00356-3  

+ Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 37, 773–782 (2019). https://doi.org/10.1038/s41587-019-0114-2 https://www.nature.com/articles/s41587-019-0114-2  

+ Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife, 6, e26476. https://doi.org/10.7554/eLife.26476  

+ Jin, H., Liu, Z. A benchmark for RNA-seq deconvolution analysis under dynamic testing environments. Genome Biol 22, 102 (2021). https://doi.org/10.1186/s13059-021-02290-6 

+ Chakravarthy, A., Furness, A., Joshi, K. et al. Pan-cancer deconvolution of tumour composition using DNA methylation. Nat Commun 9, 3220 (2018). https://doi.org/10.1038/s41467-018-05570-1 

+ Trang Le, Rachel A Aronow, Arkadz Kirshtein, Leili Shahriyari, A review of digital cytometry methods: estimating the relative abundance of cell types in a bulk of cells, Briefings in Bioinformatics, Volume 22, Issue 4, July 2021, bbaa219, https://doi.org/10.1093/bib/bbaa219 

+ Sutton, G.J., Poppe, D., Simmons, R.K. et al. Comprehensive evaluation of deconvolution methods for human brain gene expression. Nat Commun 13, 1358 (2022). https://doi.org/10.1038/s41467-022-28655-4 

+ Tran, K.A., Addala, V., Johnston, R.L. et al. Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures. Nat Commun 14, 5758 (2023). https://www.nature.com/articles/s41467-023-41385-5 

+ Vathrakokoili Pournara A, Miao Z, Beker OY, Nolte N, Brazma A, Papatheodorou I. CATD: a reproducible pipeline for selecting cell-type deconvolution methods across tissues. Bioinform Adv. 2024 Mar 23;4(1):vbae048. doi: 10.1093/bioadv/vbae048. PMID: 38638280; PMCID: PMC11023940. (EMBL)

+ Dietrich A., Metorro L., Pelz K., Eder B., Zackl C., Marini F., Sturm G., List M., Finotello F. (2024).
Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data
bioRxiv doi: https://doi.org/10.1101/2024.06.10.598226

+ Sutton, G.J., Poppe, D., Simmons, R.K. et al. Comprehensive evaluation of deconvolution methods for human brain gene expression. Nat Commun 13, 1358 (2022). https://doi.org/10.1038/s41467-022-28655-4 https://www.nature.com/articles/s41467-022-28655-4

+ Steele, N. G., Carpenter, E. S., Kemp, S. B., Sirihorachai, V. R., The, S., Delrosario, L., Lazarus, J., Amir, E. D., Gunchick, V., Espinoza, C., Bell, S., Harris, L., Lima, F., Irizarry-Negron, V., Paglia, D., Macchia, J., Chu, A. K. Y., Schofield, H., Wamsteker, E. J., Kwon, R., … Pasca di Magliano, M. (2020). Multimodal Mapping of the Tumor and Peripheral Blood Immune Landscape in Human Pancreatic Cancer. Nature cancer, 1(11), 1097–1112. https://doi.org/10.1038/s43018-020-00121-4
https://www.nature.com/articles/s43018-020-00121-4 

+ Luca, B. A., Steen, C. B., Matusiak, M., Azizi, A., Varma, S., Zhu, C., Przybyl, J., Espín-Pérez, A., Diehn, M., Alizadeh, A. A., van de Rijn, M., Gentles, A. J., & Newman, A. M. (2021). Atlas of clinically distinct cell states and ecosystems across human solid tumors. Cell, 184(21), 5482–5496.e28. https://doi.org/10.1016/j.cell.2021.09.014

+ Alexander Dietrich, Gregor Sturm, Lorenzo Merotto, Federico Marini, Francesca Finotello, Markus List, SimBu: bias-aware simulation of bulk RNA-seq data with variable cell-type composition, Bioinformatics, Volume 38, Issue Supplement_2, September 2022, Pages ii141–ii147, https://doi.org/10.1093/bioinformatics/btac499
https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii141/6702009
 
+ Benchmarking second-generation methods for cell-type deconvolution of transcriptomic data. Alexander Dietrich, Lorenzo Merotto, Konstantin Pelz, Bernhard Eder, Constantin Zackl, Katharina Reinisch, Frank Edenhofer, Federico Marini, Gregor Sturm, Markus List, Francesca Finotello. bioRxiv 2024.06.10.598226; doi: https://doi.org/10.1101/2024.06.10.598226




