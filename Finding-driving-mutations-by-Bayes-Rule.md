Finding top 10 driving mutations by Bayes’ Rule
================
2019 9 16

In this assessment, I analyzed NCBI Breast Cancer dataset by Molecular
subtypes. The description for the dataset is given below.

## DATASET

**NCBI Build 37 (UCSC hg19)**  
– Citation: Assembly \[Internet\]. Bethesda (MD): National Library of
Medicine (US), National Center for Biotechnology Information; 2012 –
\[cited 2019 09 16\]. Available from:
<https://www.ncbi.nlm.nih.gov/assembly/>

  - containing clinical information of 1,055 samples, 90,490 mutations,
    expression data of 17,814 genes for 527 samples  
  - Data was obtained with Illumina GAIIx sequencer with whole exome
    sequencing method

#### Mutation

• **Sample\_id**: mappable unique sample id of patient.  
• **Hugo\_Symbol**: gene symbol provided by HGNC  
• **ensembl\_gene\_id**: gene id provided by ensembl of EBI  
• **Chromosome**: position of each gene on chromosomes  
• **Variant\_Classification**: classification of variants.  
• **Variant\_Type**: Deletion, Insertion, or SNP  
• **Mutation\_Status**: Somatic or Germline mutation.

#### Clinical

• **sample\_id**: mappable unique sample id of patient.  
• **stage**: pathological state of each patient.  
• **subtype**: subtype of each patient classified by PAM50
classification system.  
• **survival\_time**: survival time of each patient after diagnosis of
the disease.  
• **vital\_status**: vital status of each patient. 1: died, 0: alive.  
• **survival\_index**: a class label related to the survival of the
patient and could be regarded as a class label of the same type as the
stage or subtype.

-----

## Assessment

> **1. Calculating 5-year survival probabilities of each subtypes in the
> clinical dataset.**

Load datasets.

``` r
clinical <- readRDS("./data/clinical.rds")
mutation <- readRDS("./data/mutation.rds")
```

Get the list of non-NA subtypes.

``` r
clinical_subtypes <- unique(clinical$subtype[which(!is.na(clinical$subtype))])
n_subtypes <- length(clinical_subtypes)
clinical_subtypes
```

    ## [1] Luminal A     Basal-like    Luminal B     HER2-enriched Normal-like  
    ## Levels: Basal-like HER2-enriched Luminal A Luminal B Normal-like

Calculate the 5-year survival probabilities for each subtype except
“Normal-type”. First filter out the **clinical dataset** for each
subtype, then calculate the 5-year survival probability.

``` r
for (i in 1:n_subtypes){
  if (clinical_subtypes[i] != "Normal-like"){
  clinical_subset <- dplyr::filter(clinical, subtype==clinical_subtypes[i])
  tot <- length(which(clinical_subset$survival_time>=0))
  survived <- length(which(clinical_subset$survival_time>1825))
  prob <- survived/tot
  print(paste("5-year survival prob. for ", clinical_subtypes[i], " is: ", round(prob,3)))}
}
```

    ## [1] "5-year survival prob. for  Luminal A  is:  0.301"
    ## [1] "5-year survival prob. for  Basal-like  is:  0.337"
    ## [1] "5-year survival prob. for  Luminal B  is:  0.231"
    ## [1] "5-year survival prob. for  HER2-enriched  is:  0.178"

> **2. Finding Top 10 highly possible driving mutations with estimated
> probability using conditional probability for each subtype in the
> clinical dataset.**

#### 1\. Estimation using Bayes’ Rule

For highly possible cancer-driving mutation genes of each subtype,
![p(subtype|gene)](https://latex.codecogs.com/png.latex?p%28subtype%7Cgene%29
"p(subtype|gene)") values would be high.
![p(subtype|gene)](https://latex.codecogs.com/png.latex?p%28subtype%7Cgene%29
"p(subtype|gene)") value can be calculated with Bayes’ rule.  
  
![&#10;p(S\_i|gene) = \\frac{p(gene|S\_i)p(S\_i)}{\\sum\_{i=1}^{5}
p(gene|S\_i)p(S\_i)}&#10;](https://latex.codecogs.com/png.latex?%0Ap%28S_i%7Cgene%29%20%3D%20%5Cfrac%7Bp%28gene%7CS_i%29p%28S_i%29%7D%7B%5Csum_%7Bi%3D1%7D%5E%7B5%7D%20p%28gene%7CS_i%29p%28S_i%29%7D%0A
"
p(S_i|gene) = \\frac{p(gene|S_i)p(S_i)}{\\sum_{i=1}^{5} p(gene|S_i)p(S_i)}
")  

Also, following components could be calculated as below.  
For each subtype ![S\_i](https://latex.codecogs.com/png.latex?S_i
"S_i"),

  
![&#10;p(S\_i) = \\frac{\\textbf{\# samples of subtype
i}}{\\textbf{total \#
samples}}&#10;](https://latex.codecogs.com/png.latex?%0Ap%28S_i%29%20%3D%20%5Cfrac%7B%5Ctextbf%7B%23%20samples%20of%20subtype%20i%7D%7D%7B%5Ctextbf%7Btotal%20%23%20samples%7D%7D%0A
"
p(S_i) = \\frac{\\textbf{# samples of subtype i}}{\\textbf{total # samples}}
")  

  
![&#10;p(gene\_j|S\_i) = \\frac{\\textbf{\# samples of subtype i that
have mutation j}}{\\textbf{\# samples of subtype
i}}&#10;](https://latex.codecogs.com/png.latex?%0Ap%28gene_j%7CS_i%29%20%3D%20%5Cfrac%7B%5Ctextbf%7B%23%20samples%20of%20subtype%20i%20that%20have%20mutation%20j%7D%7D%7B%5Ctextbf%7B%23%20samples%20of%20subtype%20i%7D%7D%0A
"
p(gene_j|S_i) = \\frac{\\textbf{# samples of subtype i that have mutation j}}{\\textbf{# samples of subtype i}}
")  

  
![&#10;p(S\_i|gene\_j) = \\frac{\\textbf{\# samples of subtype i that
have mutation j}{total \#
samples}&#10;](https://latex.codecogs.com/png.latex?%0Ap%28S_i%7Cgene_j%29%20%3D%20%5Cfrac%7B%5Ctextbf%7B%23%20samples%20of%20subtype%20i%20that%20have%20mutation%20j%7D%7Btotal%20%23%20samples%7D%0A
"
p(S_i|gene_j) = \\frac{\\textbf{# samples of subtype i that have mutation j}{total # samples}
")  

Find total number of sample and number of samples of each subtype.

``` r
total_sample <- unique(clinical$sample_id)
total_num_sample <- length(total_sample)
print(paste("total number of sample is ", total_num_sample))
```

    ## [1] "total number of sample is  1055"

``` r
for (i in 1:n_subtypes){
  # calculate p(S_i)
  subset <- clinical$sample_id[which(clinical$subtype==clinical_subtypes[i])]
  assign(paste("subset_",i,sep=""),dplyr::filter(mutation, sample_id %in% subset))
  p_s <- length(subset)/total_num_sample
  print(paste("p(S) for ",clinical_subtypes[i]," is: ", round(p_s,3)))}
```

    ## [1] "p(S) for  Luminal A  is:  0.186"
    ## [1] "p(S) for  Basal-like  is:  0.079"
    ## [1] "p(S) for  Luminal B  is:  0.099"
    ## [1] "p(S) for  HER2-enriched  is:  0.043"
    ## [1] "p(S) for  Normal-like  is:  0.006"

``` r
subset_na<-clinical$sample_id[which(is.na(clinical$subtype))]
assign("subset_NA",dplyr::filter(mutation, sample_id %in% subset_na))
p_s <- length(subset_na)/total_num_sample
print(paste("p(S) for nontype is: ", round(p_s,3)))
```

    ## [1] "p(S) for nontype is:  0.589"

``` r
# obtain list of mutated genes
mut_gene_list <- unique(mutation$ensembl_gene_id)
n_genes <- length(mut_gene_list)

prob_result_1 <- data.frame(gene=character(), Luminal_A=double(), Luminal_B=double(), Basal_like=double(), HER2_enriched=double(), normal_like=double(), stringsAsFactors = FALSE)

for (j in 1:n_genes){
  gene_id <- mut_gene_list[j]
  prob_result_1[j,1] <- gene_id
  #calculate sigma[p(gene|S_i)*p(S_i)]
  odd_1 <- length(which(subset_1$ensembl_gene_id==gene_id))
  odd_2 <- length(which(subset_2$ensembl_gene_id==gene_id))
  odd_3 <- length(which(subset_3$ensembl_gene_id==gene_id))
  odd_4 <- length(which(subset_4$ensembl_gene_id==gene_id))
  odd_5 <- length(which(subset_5$ensembl_gene_id==gene_id))
  odd_na <- length(which(subset_NA$ensembl_gene_id==gene_id))
  sigma <- (odd_1+odd_2+odd_3+odd_4+odd_5+odd_na)/total_num_sample
  
  # calculate p(S_i|gene) from Bayes' rule
  divisor <- total_num_sample*sigma
  prob_result_1[j,2] <- odd_1/divisor
  prob_result_1[j,3] <- odd_2/divisor
  prob_result_1[j,4] <- odd_3/divisor
  prob_result_1[j,5] <- odd_4/divisor
  prob_result_1[j,6] <- odd_5/divisor
}
```

#### 2\. Direct calculation of p(subtype|gene)

We can directly calculate
![p(S\_i|gene\_j)](https://latex.codecogs.com/png.latex?p%28S_i%7Cgene_j%29
"p(S_i|gene_j)") without using Bayes’ rule as follows.  
  
![&#10;p(S\_i|gene\_j)=\\frac{\\textbf{\# samples of subtype i with gene
j mutated}}{\\textbf{\# samples with gene j
mutated}}&#10;](https://latex.codecogs.com/png.latex?%0Ap%28S_i%7Cgene_j%29%3D%5Cfrac%7B%5Ctextbf%7B%23%20samples%20of%20subtype%20i%20with%20gene%20j%20mutated%7D%7D%7B%5Ctextbf%7B%23%20samples%20with%20gene%20j%20mutated%7D%7D%0A
"
p(S_i|gene_j)=\\frac{\\textbf{# samples of subtype i with gene j mutated}}{\\textbf{# samples with gene j mutated}}
")  

``` r
prob_result_2 <- data.frame(gene=character(), Luminal_A=double(), Luminal_B=double(), Basal_like=double(), HER2_enriched=double(), normal_like=double(), stringsAsFactors = FALSE)

for (j in 1:n_genes){
  gene_id <- mut_gene_list[j]
  # calculate the number of samples with mutation gene j
  sample_with_mut <- dplyr::filter(mutation, ensembl_gene_id==gene_id)
  sample_id_with_mut <- unique(sample_with_mut$sample_id)
  clinical_filtered<-dplyr::filter(clinical,sample_id %in% sample_id_with_mut)
  tot_sample <- nrow(clinical_filtered)
  prob_result_2[j,1] <- gene_id
  
  # calculate the number of samples with mutation gene j and subtype i 
  for (i in 1:n_subtypes){
    subtype_sample <- length(which(clinical_filtered$subtype==clinical_subtypes[i]))
    prob <- subtype_sample/tot_sample
    prob_result_2[j,1+i] <- prob
  }
}
```

Obviously, Method 1(Bayes’ rule) and Method 2(direct calculation) show
same result. However, calculation of Method 1 took less time.

#### Method 1 result

``` r
head(prob_result_1,20)
```

    ##               gene  Luminal_A  Luminal_B Basal_like HER2_enriched
    ## 1  ENSG00000014257 0.33333333 0.16666667 0.00000000    0.00000000
    ## 2  ENSG00000101901 0.00000000 0.16666667 0.16666667    0.00000000
    ## 3  ENSG00000243480 0.16666667 0.16666667 0.00000000    0.00000000
    ## 4  ENSG00000086062 0.25000000 0.25000000 0.25000000    0.00000000
    ## 5  ENSG00000132357 0.10000000 0.40000000 0.00000000    0.00000000
    ## 6  ENSG00000105479 0.40000000 0.20000000 0.00000000    0.00000000
    ## 7  ENSG00000039068 0.09917355 0.02479339 0.09917355    0.04958678
    ## 8  ENSG00000167193 0.00000000 0.00000000 0.00000000    0.00000000
    ## 9  ENSG00000215908 0.16923077 0.07692308 0.04615385    0.03076923
    ## 10 ENSG00000102034 0.12500000 0.12500000 0.12500000    0.00000000
    ## 11 ENSG00000183495 0.28571429 0.00000000 0.14285714    0.00000000
    ## 12 ENSG00000106125 0.33333333 0.00000000 0.33333333    0.11111111
    ## 13 ENSG00000107485 0.31683168 0.05940594 0.08910891    0.07920792
    ## 14 ENSG00000069122 0.23076923 0.15384615 0.07692308    0.00000000
    ## 15 ENSG00000112972 0.00000000 0.22222222 0.22222222    0.00000000
    ## 16 ENSG00000100209 0.00000000 0.00000000 0.00000000    0.00000000
    ## 17 ENSG00000124313 0.13333333 0.13333333 0.20000000    0.00000000
    ## 18 ENSG00000255103 0.18181818 0.09090909 0.18181818    0.36363636
    ## 19 ENSG00000250423 0.31034483 0.10344828 0.13793103    0.03448276
    ## 20 ENSG00000131409 0.14285714 0.28571429 0.14285714    0.00000000
    ##    normal_like
    ## 1  0.000000000
    ## 2  0.000000000
    ## 3  0.000000000
    ## 4  0.000000000
    ## 5  0.000000000
    ## 6  0.000000000
    ## 7  0.008264463
    ## 8  0.000000000
    ## 9  0.000000000
    ## 10 0.000000000
    ## 11 0.000000000
    ## 12 0.000000000
    ## 13 0.009900990
    ## 14 0.000000000
    ## 15 0.000000000
    ## 16 0.000000000
    ## 17 0.000000000
    ## 18 0.000000000
    ## 19 0.000000000
    ## 20 0.000000000

#### Method 2 result

``` r
head(prob_result_2,20)
```

    ##               gene  Luminal_A  Luminal_B Basal_like HER2_enriched
    ## 1  ENSG00000014257 0.33333333 0.16666667 0.00000000    0.00000000
    ## 2  ENSG00000101901 0.00000000 0.16666667 0.16666667    0.00000000
    ## 3  ENSG00000243480 0.16666667 0.16666667 0.00000000    0.00000000
    ## 4  ENSG00000086062 0.25000000 0.25000000 0.25000000    0.00000000
    ## 5  ENSG00000132357 0.12500000 0.25000000 0.00000000    0.00000000
    ## 6  ENSG00000105479 0.40000000 0.20000000 0.00000000    0.00000000
    ## 7  ENSG00000039068 0.09821429 0.02678571 0.08928571    0.03571429
    ## 8  ENSG00000167193 0.00000000 0.00000000 0.00000000    0.00000000
    ## 9  ENSG00000215908 0.17741935 0.08064516 0.04838710    0.03225806
    ## 10 ENSG00000102034 0.14285714 0.14285714 0.14285714    0.00000000
    ## 11 ENSG00000183495 0.17647059 0.00000000 0.11764706    0.00000000
    ## 12 ENSG00000106125 0.28571429 0.00000000 0.28571429    0.14285714
    ## 13 ENSG00000107485 0.31914894 0.05319149 0.09574468    0.06382979
    ## 14 ENSG00000069122 0.18181818 0.09090909 0.09090909    0.00000000
    ## 15 ENSG00000112972 0.00000000 0.12500000 0.25000000    0.00000000
    ## 16 ENSG00000100209 0.00000000 0.00000000 0.00000000    0.00000000
    ## 17 ENSG00000124313 0.14285714 0.14285714 0.21428571    0.00000000
    ## 18 ENSG00000255103 0.25000000 0.12500000 0.25000000    0.12500000
    ## 19 ENSG00000250423 0.19047619 0.04761905 0.19047619    0.04761905
    ## 20 ENSG00000131409 0.14285714 0.28571429 0.14285714    0.00000000
    ##    normal_like
    ## 1  0.000000000
    ## 2  0.000000000
    ## 3  0.000000000
    ## 4  0.000000000
    ## 5  0.000000000
    ## 6  0.000000000
    ## 7  0.008928571
    ## 8  0.000000000
    ## 9  0.000000000
    ## 10 0.000000000
    ## 11 0.000000000
    ## 12 0.000000000
    ## 13 0.010638298
    ## 14 0.000000000
    ## 15 0.000000000
    ## 16 0.000000000
    ## 17 0.000000000
    ## 18 0.000000000
    ## 19 0.000000000
    ## 20 0.000000000

## Result

#### TOP 10 driving mutations for each subtype

``` r
dplyr::select(head(prob_result_1[order(-prob_result_1$Luminal_A),],10),gene,Luminal_A)
```

    ##                gene Luminal_A
    ## 406 ENSG00000136932         1
    ## 410 ENSG00000213088         1
    ## 444 ENSG00000167333         1
    ## 451 ENSG00000064666         1
    ## 452 ENSG00000105219         1
    ## 711 ENSG00000215394         1
    ## 722 ENSG00000199704         1
    ## 762 ENSG00000211648         1
    ## 833 ENSG00000197445         1
    ## 850 ENSG00000110218         1

``` r
dplyr::select(head(prob_result_1[order(-prob_result_1$Luminal_B),],10),gene,Luminal_B)
```

    ##                 gene Luminal_B
    ## 313  ENSG00000073598         1
    ## 333  ENSG00000134330         1
    ## 349  ENSG00000167965         1
    ## 389  ENSG00000196800         1
    ## 1030 ENSG00000105173         1
    ## 1889 ENSG00000186453         1
    ## 1994 ENSG00000171720         1
    ## 1996 ENSG00000211944         1
    ## 2015 ENSG00000175792         1
    ## 3208 ENSG00000158525         1

``` r
dplyr::select(head(prob_result_1[order(-prob_result_1$Basal_like),],10),gene,Basal_like)
```

    ##                 gene Basal_like
    ## 1894 ENSG00000103365          1
    ## 1935 ENSG00000136868          1
    ## 1953 ENSG00000184560          1
    ## 1960 ENSG00000053254          1
    ## 1962 ENSG00000174564          1
    ## 1965 ENSG00000167749          1
    ## 1978 ENSG00000116209          1
    ## 3052 ENSG00000166947          1
    ## 3065 ENSG00000140939          1
    ## 3067 ENSG00000176269          1

``` r
dplyr::select(head(prob_result_1[order(-prob_result_1$HER2_enriched),],10),gene,HER2_enriched)
```

    ##                 gene HER2_enriched
    ## 1828 ENSG00000092531             1
    ## 1941 ENSG00000253797             1
    ## 3834 ENSG00000244731             1
    ## 3835 ENSG00000204839             1
    ## 3841 ENSG00000086619             1
    ## 8345 ENSG00000106009             1
    ## 8355 ENSG00000166896             1
    ## 8420 ENSG00000182575             1
    ## 8551 ENSG00000182050             1
    ## 8737 ENSG00000169717             1

``` r
dplyr::select(head(prob_result_1[order(-prob_result_1$normal_like),],10),gene,normal_like)
```

    ##                 gene normal_like
    ## 9373 ENSG00000137090           1
    ## 9374 ENSG00000183470           1
    ## 9375 ENSG00000210194           1
    ## 9376 ENSG00000265489           1
    ## 9384 ENSG00000211961           1
    ## 9385 ENSG00000178502           1
    ## 9390 ENSG00000090971           1
    ## 9396 ENSG00000182676           1
    ## 9398 ENSG00000137522           1
    ## 9405 ENSG00000196860           1

## Discussion

Since the number of samples that has each mutated gene was very small
(1\~7), there was a high probability for all samples for each gene
having the same subtype. This made the probabilities
![p(subtype\_i|gene\_j)](https://latex.codecogs.com/png.latex?p%28subtype_i%7Cgene_j%29
"p(subtype_i|gene_j)") to be 1 so easily. Because of this limitation, it
is difficult to conclude that the shown Top 10 genes are the highest
cancer-driving mutation genes.

Method 1 took less time compared to Method 2. It would be interesting to
further compare the computational cost of R function such as
‘dplyr::filter’, ‘which’, ‘unique’, etc.

## References

  - Tutorial Manual (KAIST BiS335)
