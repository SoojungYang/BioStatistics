Support Vector Classifier with 3 kernels
================
Soojung Yang
2019 12 8

The goal of this assessment is to implement Support Vector Classifier
with various kernel functions (linear, polynomial, radial).

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: lattice

    ## Loading required package: ggplot2

#### **DATASET**

**NCBI Build 37 (UCSC hg19)**  
– Citation: Assembly \[Internet\]. Bethesda (MD): National Library of
Medicine (US), National Center for Biotechnology Information; 2012 –
\[cited 2019 09 16\]. Available from:
<https://www.ncbi.nlm.nih.gov/assembly/>

-----

#### **Data Input**

``` r
# Load data
clin <- readRDS("./data/clinical.rds")
gex <- readRDS("./data/expression.rds")
```

#### **Data preprocessing**

First, we select the clinical data points only if they are in expression
dataset. Then, we remove genes that have NA expressions. We normalize
the expression values between samples, and balance the stages by
undersampling. Finally, by performing ANOVA, we select top 10 genes for
efficient analysis.

``` r
# Clinical data preprocessing
idx <- intersect(colnames(gex), clin$sample)
clin_p <- clin[clin$sample_id %in% idx, ] # preprocessed clin
stage <- as.vector(clin_p$stage)
stage_p <- rep(0,length(stage))
for (i in 1:length(stage)){
if(stage[i]=='stage i'){stage_p[i]=1}
else if(stage[i]=='stage ii'){stage_p[i]=2}
else if(stage[i]=='stage iii'){stage_p[i]=3}
else if(stage[i]=='stage iv'){stage_p[i]=4}
else if(stage[i]=='stage v'){stage_p[i]=5}
}
## Gene expression data preprocessing
gex_p <- gex[, colnames(gex) %in% idx]
gex_p <- na.omit(gex_p) # Remove genes containing NA

# Quantile Normalize (Between sample normalization)
temp <- normalize.quantiles(gex_p)
colnames(temp) <- colnames(gex_p)
rownames(temp) <- rownames(gex_p)
gex_p <- temp
# Preprocessed data
data <- data.frame(stage_p, t(gex_p))
# Preprocessing for solving class imbalance (Undersampling)
set.seed(1)
stage2.sam <- sample(rownames(data)[data$stage_p == 2], size=min(sum(data$stage_p==1)))
stage3.sam <- sample(rownames(data)[data$stage_p == 3], size=min(sum(data$stage_p==1)))
data <- data[rownames(data) %in% c(rownames(data)[data$stage_p == 1], stage2.sam, stage3.sam), ]
# Gene selection for decreasing computation time using ANOVA
feature.size <- 10
aov.gene.pval <- c()
for(i in 2:(ncol(data))){
aov.gene.pval <- c(aov.gene.pval, summary(aov(data$stage_p ~ data[, i]))[[1]][["Pr(>F)"]][1])
}
aov.idx <- order(aov.gene.pval, decreasing=FALSE)[1:feature.size] # (ANOVA P value top n genes)
data.aov <- data[, c(1, 1+aov.idx)] # 1st column is stage_p
```

Dataset partioning was done by createDataPartition function, with fixed
random seed.

``` r
set.seed(1)
# stage_p numeric to factor for classification
data.aov$stage_p <- as.factor(data.aov$stage_p)
# training:test = 0.8:0.2
partition <- createDataPartition(data.aov$stage_p, p=0.8, list=FALSE)
train <- data.aov[partition,]
test <- data.aov[-partition,]
```

### Linear kernel SVM

``` r
linear.tune <- tune(svm, stage_p~., data = train, kernel = "linear", ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10)))
summary(linear.tune)
```

    ## 
    ## Parameter tuning of 'svm':
    ## 
    ## - sampling method: 10-fold cross validation 
    ## 
    ## - best parameters:
    ##  cost
    ##   0.1
    ## 
    ## - best performance: 0.5512987 
    ## 
    ## - Detailed performance results:
    ##    cost     error dispersion
    ## 1 1e-03 0.7456710 0.11365526
    ## 2 1e-02 0.5878788 0.08766278
    ## 3 1e-01 0.5512987 0.08667653
    ## 4 1e+00 0.5748918 0.09152088
    ## 5 1e+01 0.5839827 0.08476820

``` r
linear.model <- linear.tune$best.model
summary(linear.model)
```

    ## 
    ## Call:
    ## best.tune(method = svm, train.x = stage_p ~ ., data = train, 
    ##     ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10)), kernel = "linear")
    ## 
    ## 
    ## Parameters:
    ##    SVM-Type:  C-classification 
    ##  SVM-Kernel:  linear 
    ##        cost:  0.1 
    ## 
    ## Number of Support Vectors:  199
    ## 
    ##  ( 64 72 63 )
    ## 
    ## 
    ## Number of Classes:  3 
    ## 
    ## Levels: 
    ##  1 2 3

``` r
linear.pred <- predict(linear.model, test)
confusionMatrix(linear.pred, test$stage_p)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 11  3  3
    ##          2  6  4  5
    ##          3  0 10  9
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.4706          
    ##                  95% CI : (0.3293, 0.6154)
    ##     No Information Rate : 0.3333          
    ##     P-Value [Acc > NIR] : 0.02892         
    ##                                           
    ##                   Kappa : 0.2059          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.12900         
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.6471  0.23529   0.5294
    ## Specificity            0.8235  0.67647   0.7059
    ## Pos Pred Value         0.6471  0.26667   0.4737
    ## Neg Pred Value         0.8235  0.63889   0.7500
    ## Prevalence             0.3333  0.33333   0.3333
    ## Detection Rate         0.2157  0.07843   0.1765
    ## Detection Prevalence   0.3333  0.29412   0.3725
    ## Balanced Accuracy      0.7353  0.45588   0.6176

The accuracy for linear kernel SVM model is

``` r
mean(linear.pred == test$stage_p)
```

    ## [1] 0.4705882

### Polynomial kernel SVM

``` r
# Polynomial kernel SVM 
poly.tune <- tune(svm, stage_p~., data = train, kernel = "polynomial", ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10), gamma = c(0.01, 0.5, 1, 2)))
summary(poly.tune)
```

    ## 
    ## Parameter tuning of 'svm':
    ## 
    ## - sampling method: 10-fold cross validation 
    ## 
    ## - best parameters:
    ##  cost gamma
    ##  0.01   0.5
    ## 
    ## - best performance: 0.6015152 
    ## 
    ## - Detailed performance results:
    ##     cost gamma     error dispersion
    ## 1  1e-03  0.01 0.7406926 0.04944044
    ## 2  1e-02  0.01 0.7406926 0.04944044
    ## 3  1e-01  0.01 0.7406926 0.04944044
    ## 4  1e+00  0.01 0.7406926 0.04944044
    ## 5  1e+01  0.01 0.7220779 0.06600214
    ## 6  1e-03  0.50 0.6800866 0.08242052
    ## 7  1e-02  0.50 0.6015152 0.08069495
    ## 8  1e-01  0.50 0.6755411 0.10737254
    ## 9  1e+00  0.50 0.6852814 0.06703514
    ## 10 1e+01  0.50 0.6673160 0.07077367
    ## 11 1e-03  1.00 0.6246753 0.08074492
    ## 12 1e-02  1.00 0.6621212 0.12019626
    ## 13 1e-01  1.00 0.6714286 0.05849765
    ## 14 1e+00  1.00 0.6675325 0.08248618
    ## 15 1e+01  1.00 0.6673160 0.06711935
    ## 16 1e-03  2.00 0.6571429 0.11325125
    ## 17 1e-02  2.00 0.6987013 0.05259339
    ## 18 1e-01  2.00 0.6629870 0.07348016
    ## 19 1e+00  2.00 0.6673160 0.06711935
    ## 20 1e+01  2.00 0.6673160 0.06711935

``` r
poly.model <- poly.tune$best.model
summary(poly.model)
```

    ## 
    ## Call:
    ## best.tune(method = svm, train.x = stage_p ~ ., data = train, 
    ##     ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10), gamma = c(0.01, 
    ##         0.5, 1, 2)), kernel = "polynomial")
    ## 
    ## 
    ## Parameters:
    ##    SVM-Type:  C-classification 
    ##  SVM-Kernel:  polynomial 
    ##        cost:  0.01 
    ##      degree:  3 
    ##      coef.0:  0 
    ## 
    ## Number of Support Vectors:  210
    ## 
    ##  ( 67 72 71 )
    ## 
    ## 
    ## Number of Classes:  3 
    ## 
    ## Levels: 
    ##  1 2 3

``` r
poly.pred <- predict(poly.model, test)
confusionMatrix(poly.pred, test$stage_p)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 12  2  4
    ##          2  4  4  3
    ##          3  1 11 10
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.5098         
    ##                  95% CI : (0.366, 0.6525)
    ##     No Information Rate : 0.3333         
    ##     P-Value [Acc > NIR] : 0.006888       
    ##                                          
    ##                   Kappa : 0.2647         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.070693       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.7059  0.23529   0.5882
    ## Specificity            0.8235  0.79412   0.6471
    ## Pos Pred Value         0.6667  0.36364   0.4545
    ## Neg Pred Value         0.8485  0.67500   0.7586
    ## Prevalence             0.3333  0.33333   0.3333
    ## Detection Rate         0.2353  0.07843   0.1961
    ## Detection Prevalence   0.3529  0.21569   0.4314
    ## Balanced Accuracy      0.7647  0.51471   0.6176

The accuracy for polynomial kernel SVM model is

``` r
mean(poly.pred == test$stage_p)
```

    ## [1] 0.5098039

### Radial kernel SVM

``` r
# Radial kernel SVM
rad.tune <- tune(svm, stage_p~., data = train, kernel = "radial", ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10, 100), gamma = c(0.01, 0.5, 1, 2, 4)))
summary(rad.tune)
```

    ## 
    ## Parameter tuning of 'svm':
    ## 
    ## - sampling method: 10-fold cross validation 
    ## 
    ## - best parameters:
    ##  cost gamma
    ##    10  0.01
    ## 
    ## - best performance: 0.5603896 
    ## 
    ## - Detailed performance results:
    ##     cost gamma     error dispersion
    ## 1  1e-03  0.01 0.7774892 0.05399988
    ## 2  1e-02  0.01 0.7774892 0.05399988
    ## 3  1e-01  0.01 0.7590909 0.07530650
    ## 4  1e+00  0.01 0.5746753 0.10009506
    ## 5  1e+01  0.01 0.5603896 0.06359406
    ## 6  1e+02  0.01 0.6155844 0.07353788
    ## 7  1e-03  0.50 0.7774892 0.05399988
    ## 8  1e-02  0.50 0.7774892 0.05399988
    ## 9  1e-01  0.50 0.7774892 0.05399988
    ## 10 1e+00  0.50 0.6560606 0.13307904
    ## 11 1e+01  0.50 0.6551948 0.13544646
    ## 12 1e+02  0.50 0.6551948 0.13544646
    ## 13 1e-03  1.00 0.7774892 0.05399988
    ## 14 1e-02  1.00 0.7774892 0.05399988
    ## 15 1e-01  1.00 0.7774892 0.05399988
    ## 16 1e+00  1.00 0.6939394 0.11000082
    ## 17 1e+01  1.00 0.6930736 0.11499191
    ## 18 1e+02  1.00 0.6930736 0.11499191
    ## 19 1e-03  2.00 0.7774892 0.05399988
    ## 20 1e-02  2.00 0.7774892 0.05399988
    ## 21 1e-01  2.00 0.7774892 0.05399988
    ## 22 1e+00  2.00 0.7502165 0.09160616
    ## 23 1e+01  2.00 0.7318182 0.09569858
    ## 24 1e+02  2.00 0.7318182 0.09569858
    ## 25 1e-03  4.00 0.7774892 0.05399988
    ## 26 1e-02  4.00 0.7774892 0.05399988
    ## 27 1e-01  4.00 0.7774892 0.05399988
    ## 28 1e+00  4.00 0.7636364 0.06464421
    ## 29 1e+01  4.00 0.7545455 0.06943821
    ## 30 1e+02  4.00 0.7545455 0.06943821

``` r
rad.model <- rad.tune$best.model
summary(rad.model)
```

    ## 
    ## Call:
    ## best.tune(method = svm, train.x = stage_p ~ ., data = train, 
    ##     ranges = list(cost = c(0.001, 0.01, 0.1, 1, 10, 100), gamma = c(0.01, 
    ##         0.5, 1, 2, 4)), kernel = "radial")
    ## 
    ## 
    ## Parameters:
    ##    SVM-Type:  C-classification 
    ##  SVM-Kernel:  radial 
    ##        cost:  10 
    ## 
    ## Number of Support Vectors:  200
    ## 
    ##  ( 66 72 62 )
    ## 
    ## 
    ## Number of Classes:  3 
    ## 
    ## Levels: 
    ##  1 2 3

``` r
rad.pred <- predict(rad.model, test)
confusionMatrix(rad.pred, test$stage_p)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3
    ##          1 11  1  4
    ##          2  6  7  6
    ##          3  0  9  7
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.4902         
    ##                  95% CI : (0.3475, 0.634)
    ##     No Information Rate : 0.3333         
    ##     P-Value [Acc > NIR] : 0.01461        
    ##                                          
    ##                   Kappa : 0.2353         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.04260        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3
    ## Sensitivity            0.6471   0.4118   0.4118
    ## Specificity            0.8529   0.6471   0.7353
    ## Pos Pred Value         0.6875   0.3684   0.4375
    ## Neg Pred Value         0.8286   0.6875   0.7143
    ## Prevalence             0.3333   0.3333   0.3333
    ## Detection Rate         0.2157   0.1373   0.1373
    ## Detection Prevalence   0.3137   0.3725   0.3137
    ## Balanced Accuracy      0.7500   0.5294   0.5735

The accuracy for radial kernel SVM model is

``` r
mean(rad.pred == test$stage_p)
```

    ## [1] 0.4901961

## **Assessment**

> 1.  Build svm model changing kernel function to classify patient
>     stage. And compare their performance.

The accuracy of the model were 0.4706 for linear kernel, 0.5098 for
polynomial kernel, and 0.4902 for radial kernel. Comparing the accuracy,
polynomial kernel model showed the best performance, and linear kernel
model was the worst. Other metrics such as sensitivity, specificity,
Rand index, or Dice index can be used to compare the model. In this
assessment, I compared the models with simple accuracy since the test
set is well-balanced between classes(stage). Also, we can qualitatively
compare the model with confusion matrix. The reason why the linear
kernel showed the worst performance might be that it was too simple for
this data.

> 2.  In case of radial kernel model, find its optimal cost & gamma
>     parameter. Wrtie its mathematical description and relation with
>     model’s bias-variance.

The optimization object for SVM is  
![min \\frac12w^Tw +
C\\sum\\limits\_{i=1}^{p}e\_i](https://latex.codecogs.com/png.latex?min%20%5Cfrac12w%5ETw%20%2B%20C%5Csum%5Climits_%7Bi%3D1%7D%5E%7Bp%7De_i
"min \\frac12w^Tw + C\\sum\\limits_{i=1}^{p}e_i") Subject to
![y\_i(w^T\\phi(x\_i)+b)\\ge1-e\_i,
e\_i\\ge0](https://latex.codecogs.com/png.latex?y_i%28w%5ET%5Cphi%28x_i%29%2Bb%29%5Cge1-e_i%2C%20e_i%5Cge0
"y_i(w^T\\phi(x_i)+b)\\ge1-e_i, e_i\\ge0")  
Also, for radial kernel function, ![K(x\_i,
x\_j)](https://latex.codecogs.com/png.latex?K%28x_i%2C%20x_j%29
"K(x_i, x_j)") for
![\\phi(x)](https://latex.codecogs.com/png.latex?%5Cphi%28x%29
"\\phi(x)") is ![K(x\_i, x\_j) = exp(-\\gamma\\parallel x\_i - x\_j
\\parallel^2), \\gamma
\> 0](https://latex.codecogs.com/png.latex?K%28x_i%2C%20x_j%29%20%3D%20exp%28-%5Cgamma%5Cparallel%20x_i%20-%20x_j%20%5Cparallel%5E2%29%2C%20%5Cgamma%20%3E%200
"K(x_i, x_j) = exp(-\\gamma\\parallel x_i - x_j \\parallel^2), \\gamma \> 0").  
The cost is C from the first inequality, and gamma is
![\\gamma](https://latex.codecogs.com/png.latex?%5Cgamma "\\gamma") from
the second inequality.  
The optimal cost value was 10 and gamma value was 0.01.  
Cost(C) is the penalty parameter of the error term. Lower the cost,
larger the possible
margin(![e\_i](https://latex.codecogs.com/png.latex?e_i "e_i")), and
higher the bias & lower the variance. Larger margin means that we allow
many exceptions. Thus, the classifier becomes less fit to training data,
and the variance becomes lower. The lower variance trades off with
higher bias.  
Gamma parameter determines how far the influence from each data point
can reach the other point. From the second inequality, small gamma
implies that the class of support vector
![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") will have
influence on class of vector
![x\_j](https://latex.codecogs.com/png.latex?x_j "x_j") even when those
two points are far apart. In other words, small gamma leads to large
radius of influence of samples, and the model will be constrained to the
training data. This means the model has low bias, and high variance as a
tradeoff. By similar rational, large gamma leads to high bias and low
variance.  
If we let ![\\gamma =
\\frac{1}{2\\sigma^2}](https://latex.codecogs.com/png.latex?%5Cgamma%20%3D%20%5Cfrac%7B1%7D%7B2%5Csigma%5E2%7D
"\\gamma = \\frac{1}{2\\sigma^2}"), the radial kernel function becomes
Gaussian, or RBF(radial basis function). In this case,
![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma "\\sigma") is a
new free parameter instead of
![\\gamma](https://latex.codecogs.com/png.latex?%5Cgamma "\\gamma").
This representation gives more intuitive explanation on the role of
![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma "\\sigma"). We
can consider ![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma
"\\sigma") value as variance, so that large
![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma "\\sigma")
(which is same as small
![\\gamma](https://latex.codecogs.com/png.latex?%5Cgamma "\\gamma"))
means Gaussian with large variance, implying that the impact of each
support vector is larger.

## References

  - \[1\] Tutorial Manual (KAIST BiS335)
