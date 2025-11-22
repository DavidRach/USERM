# Basic instruction for USERM use

Xiangming Cai

2025-10-20

## 🔍 Introduction

The USERM package provides an out-of-box tool to apply the residual model approach [^1], which characterizes and predicts the spread of unmixed spectral flow cytometry data, which arises from instrumental noise or deviations between actual cellular emission and the average fluorescence signatures.

The USERM also supports computing various matrixes tools for panel design, including the Coef Matrix, the Hotspot Matrix, and the Similarity Matrix, and the Spread Distance Matrix.


## 💻 Installation

You can install the development version of USERM from [GitHub](https://github.com/xiangmingcai) with:

``` r
devtools::install_github("xiangmingcai/USERM")
```

## Step 1 📈 query avalilable fluorescence list

You can use querySig function to query all available fluorescence.

``` r
library(USERM)
Sig_info = querySig()
head(Sig_info)
```

``` r
> head(Sig_info)
                    id PrimaryName SecondaryName detectors instrument       Source Note
1  SCC_Bead_CD4_NFR700         CD4        NFR700        51     Xenith XiangmingCai   NA
2   SCC_Bead_CD3_BV510         CD3         BV510        51     Xenith XiangmingCai   NA
3    SCC_Bead_CD2_FITC         CD2          FITC        51     Xenith XiangmingCai   NA
4       SCC_Bead_AF_AF          AF            AF        51     Xenith XiangmingCai   NA
5   SCC_Bead_CD8_BV570         CD8         BV570        51     Xenith XiangmingCai   NA
6 SCC_Bead_CD16_BUV805        CD16        BUV805        51     Xenith XiangmingCai   NA
```

To check signature of a specific fluroescence, you may

``` r
fluor_to_check = Sig_info$id[3]
print(fluor_to_check)
checkSig_linePlot(id = fluor_to_check)
```
``` r
> print(fluor_to_check)
[1] "SCC_Bead_CD2_FITC"
```

<p align="center">

<img src="./images/checkSig_linePlot.jpg" width="500/"/>

</p>

These available fluorescence were preprocessed with parameters extracted for redisual model. To access the corresponding residaul model and its parameters, you can use the getRes function.

``` r
ResObj = getRes(id = fluor_to_check)
```
<p align="center">

<img src="./images/ResObj.jpg" width="500/"/>

</p>

The ResObj Object contains all relavant informations. Please read our paper (mentioned on the top) to know the concepts of each matrix shown here. See \code{\link{CreateRes}} and \code{\link{SlopEstimation}} for detailed introduction of the structure of ResObj.

You can easily visaulize all these matrixes, when you are interested in any one of them.

To check the slop matrix:
``` r
checkRes_slopMtx(Res = ResObj)
```
<p align="center">

<img src="./images/slopmtx.jpg" width="500/"/>

</p>

To check the intercept matrix:
``` r
checkRes_interceptMtx(Res = ResObj)
```
<p align="center">

<img src="./images/interceptmtx.jpg" width="500/"/>

</p>

You can check the covariance matrix at a specified fluorescence intensity (bin):
``` r
print(ResObj$bin_mids)
checkRes_covMtx(Res = ResObj,bin=3)
```
<p align="center">

<img src="./images/covmtx.jpg" width="500/"/>

</p>

You can even check the linear fitting of covariance versus fluorescence intensity (bin) between two speific detectors:
``` r
print(ResObj$detectors)
checkRes_covScatter(Res = ResObj,
                    detector1 = ResObj$detectors[1],
                    detector2 = ResObj$detectors[2])
```
<p align="center">

<img src="./images/covscatter.jpg" width="500/"/>

</p>


## Step 2 🖊 Select fluors and create SuprObj

You may now select fluorescence to make a panel. Of note, the targets (e.g. CD3, CD4) do not matter here, as long as the fluorescence matches. That means if you plan to use FITC-CD4 in you panel, you can use FITC-CD2 in the USERM to make estimation and prediction. 
If the fluorescence that you need were not available in the USERM, you can prepare single color control FCS file yourself and process it following other [instruction](../instruction_customFCS.md). 

Assume all fluorescence that you need is available, now you can select them out:
``` r
fluors_selected = c(Sig_info$id[c(32:38,41:42,46:48,50:54,8,9,12,63)])
print(fluors_selected)
```
``` r
> print(fluors_selected)
 [1] "SCC_Cell_LD_LDNIR876"       "SCC_Cell_CD2_FITC"          "SCC_Cell_CD3_BV510"        
 [4] "SCC_Cell_CD4_NFR700"        "SCC_Cell_CD8_BV570"         "SCC_Cell_CD16_BUV805"      
 [7] "SCC_Cell_CD19_BUV496"       "SCC_Cell_CD38_BV650"        "SCC_Cell_CD45_AF532"       
[10] "SCC_Cell_CD127_BUV661"      "SCC_Cell_CD161_BUV615"      "SCC_Cell_CD244_SB600"      
[13] "SCC_Cell_CRTH2_PECy5"       "SCC_Cell_HLADR_NFB61070S"   "SCC_Cell_KIRDL1_PE"        
[16] "SCC_Cell_KLRG1_BV480"       "SCC_Cell_NKG2A_PEDazzle594" "SCC_Bead_CD25_PECy55"      
[19] "SCC_Bead_CD34_APCeFluor780" "SCC_Bead_CD56_BUV563"       "SCC_Cell_AF_AF"   
```
Here we select 20 fluorescence plus autofluorescence (AF) to make a panel. Notely, both SCC Cells and Beads are using together. Based on our research, it is ok to use SCC beads for prediction.


Now, we can create the signature matrix of selected fluorescence
``` r
Sig_mtx  = getSigMtx(ids = fluors_selected) #(detectors x fluorescence)
dim(Sig_mtx)
```
``` r
> dim(Sig_mtx)
[1] 51 21
```

Then, we create a UsermObj to contain all ResObj. 

``` r
UsermObj = CreateUserm(A = Sig_mtx)
#add ResObj into UsermObj
for (save_suf in colnames(Sig_mtx)) {
  ResObj = getRes(id = save_suf)
  UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
}
```
<p align="center">

<img src="./images/UsermObj.jpg" width="500/"/>

</p>

The UsermObj contains all information needed for prediction.  The Res item is a list with ResObj of each fluorescence. The Scale_df item contains the scale setting for visualization afterward. It is ok to leave it as NA. Default linear scale will be applied (Options: "Linear", "Log10", "Arcsinh"). The Intensity_mtx contains the estimated fluorescence intensities of cell populations. In the example shown here, we have one default populaiton named "V1". You can modify the name or intensities or add new populations for interactive prediction in Step 3.

## Step 3 🧮 Make prediction (optional)

We can now interactively predict spread from residual model, which accounts for instrumental noise and deviations between actual cellular emission and the average fluorescence signatures.

To predict for one population:
``` r
PredOneSpread(Userm = UsermObj,population_id = c("V1"))
```
<p align="center">

<img src="./images/prediction_1.jpg" width="500/"/>

</p>

You will see a pop-out shinyapp window. 

<p align="center">

<img src="./images/prediction_2.jpg" width="500/"/>

</p>

The left panel allows users to interactively select and un-select fluroescence and detectors, and modify intensities. So that users can flexibily explore the impact of fluorescence combination and intensities on the predicted spread.

<p align="center">

<img src="./images/prediction_3.jpg" width="500/"/>

</p>

The right panel have some tabs. The unmixing matrix is signature matrix.

<p align="center">

<img src="./images/prediction_4.jpg" width="500/"/>

</p>

In the prediction tab, you will see a scatter plot visualizing the 95% range of predicted spread. Two display mode is supported: ***Pseudo-color*** and ***Counter-line***. It allows the selection of fluorescence for x and y axes. You can also adjsut the display limis of the axes. Scale is support with ***Linear***, ***Log10*** and ***Arcsinh*** options. For Arcsinh, users need to set a cofactor.

<p align="center">

<img src="./images/af_factor.jpg" width="500/"/>

</p>

Of note, you can now set the ***Factor for default Autofluorescence*** to control the displayed spread from AF. the displayed spread will be multiplied with the (factor / 100). So, when the factor is set to 100, the estimated spread from AF will be displayed completely. If the factor is set to 0, the spread from AF will be removed. This helps to check only the spread from noise. 

<p align="center">

<img src="./images/prediction_5.jpg" width="500/"/>

</p>

The pinv matrix shows the pseudo-inverse of the signature matrix.

<p align="center">

<img src="./images/prediction_6.jpg" width="500/"/>

</p>

In the intercept matrix tab, you can check the intercept matries of each fluroescence. 

<p align="center">

<img src="./images/prediction_7.jpg" width="500/"/>

</p>

The weighted intercept matrix for each fluorescence is calculated as <img src="./images/function1.jpg" width="100/"/> The median value of all weighted intercept matrixes is used in the prediction (last two columns).

<p align="center">

<img src="./images/prediction_8.jpg" width="500/"/>

</p>

In the slop matrix tab, you can check the slop matries of each fluroescence. 

<p align="center">

<img src="./images/prediction_9.jpg" width="500/"/>

</p>

Similarly, The weighted slop matrix for each fluorescence is calculated as <img src="./images/function2.jpg" width="100/"/> The summed value of all weighted slop matrixes is used in the prediction (last two columns).

<p align="center">

<img src="./images/prediction_10.jpg" width="500/"/>

</p>

In the export tab, you may directlyh download a complete report of all the matrixes.
You may also export the N-by-N or N-by-1 plots, which could be quite big for a big panel. If you want to plot the N-by-N or N-by-1 plots, you need to assign corresponding values to the UsermObj[["Scale_df"]]. Otherwise, an error will occur.

Prediction for more than one population is possible. However, users need to set the intensities of populations in advance. Here we predict 3 populations together.
``` r
UsermObj$Intensity_mtx[,2] = 200*c(1:length(UsermObj$fluors))
UsermObj$Intensity_mtx[,3] = 300*c(1:length(UsermObj$fluors))
names(UsermObj$Intensity_mtx) = c("V1","V2","V3")
PredMultipleSpread(Userm = UsermObj,population_ids = c("V1","V2","V3"))
```
<p align="center">

<img src="./images/prediction_11.jpg" width="500/"/>

</p>
All 3 populations will be shown together in the scatter plot. The design od the interactive web is similar to that for PredOneSpread function.

## Step 4 matrics estimation and visualization

The USERM supports multiple matrix tools for panel design, including the Coef Matrix, the Hotspot Matrix, and the Similarity Matrix, and the Spread Distance Matrix. All matrix can be visualized with the ***Vis_Mtx*** function.

To estimate a Coef Matrix:
``` r
Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Coef",
        title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")

```
<p align="center">

<img src="./images/matrix_1.jpg" width="500/"/>

</p>

In the Coef Matrix, a high coef number indicates high increasing trend of spread from one fluorescence (row) into another fluorescence (column), following the increase of fluorescence intensity. To interprate a high observed coef number and its related severe spread, researchers can use PredOneSpread function to look for the underlying machanical explanation and potentail solutions. A detailed example is provided in another [instruction](../instruction_interprateHighCoef.md). 

There is no specific threshold for "bad" coef number or "good" coef number. Because the fluorescence intensity also matters and it varies in different settings and conjugated targets. However, this Coef Matrix provides an approach to compare the estimated spread between panels. So that users can choose a panel with less expected spread. 


To estimate a Hotspot Matrix, only the signature matrix A is required:
``` r
Hotspot_mtx = EstimateHotspotMtx(A = UsermObj$A)
Vis_Mtx(mat = Hotspot_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Hotspot",
        title = "Hotspot matrix")

```
<p align="center">

<img src="./images/matrix_2.jpg" width="500/"/>

</p>

Based on the paper from Mage PL and et al.[^2], The diagonal entries of the Hotspot matrix are informative, as they reflect the degree to which each fluorescence contributes to unmixing-dependent spread.  

To estimate a Similarity Matrix, only the signature matrix A is required:
``` r
Similarity_mtx = EstimateSimilarityMtx(A = UsermObj$A)
Vis_Mtx(mat = Similarity_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 1,mid = 0.8,min = 0,legend_name = "Cosine",
        title = "Cosine similarity matrix")
```
<p align="center">

<img src="./images/matrix_3.jpg" width="500/"/>

</p>

For similarity matrix, a number higher than 0.98 indicates that signals from the pair of fluorescence cannot be well separated. So, one of them needs to be replaced by other fluorescence.

It is also possible to estimate the Spread Distance Matrix between two populations with estimated fluorescence intensities:
``` r
Spr1 = EstimateSpread(Userm = UsermObj,population_id = c("P1"))
Spr2 = EstimateSpread(Userm = UsermObj,population_id = c("P2"))
SpreadDistance_mtx = EstimateDistance(Spr1, Spr2)
Vis_Mtx(mat = SpreadDistance_mtx,mincolor = "darkred",midcolor = "white", maxcolor = "white",max = 1.5,mid = 1,min = 0,title = "Spread Distance matrix")
```
<p align="center">

<img src="./images/matrix_4.jpg" width="500/"/>

</p>

For the Spread Distance matrix, a value < 1.0 indicates low resolution between these two cell populations. Of note, this results only accounts for spread, which originates from instrumental noise or deviations between actual cellular emission and the average fluorescence signatures. Also, it only holds for the given fluorescence intensities.

## 📚 Citation

If you use this package in your research, please cite our paper and the package as:
```
Xiangming Cai, Sara Garcia-Garcia, Michaela Gianniou, Juan J. Garcia Vallejo. Manuscript in preparation. (to be update)

Cai X (2025). _USERM: Unmixing Spread Estimation with Residual Model_. R package version
  0.0.0.9000,  <https://github.com/xiangmingcai/USERM>.
  
@Manual{,
    title = {USERM: Unmixing Spread Estimation with Residual Model},
    author = {Xiangming Cai},
    year = {2025},
    note = {R package version 0.1},
    url = {https://github.com/xiangmingcai/USERM},
  }
```




[^1]: Xiangming Cai, Sara Garcia-Garcia, Michaela Gianniou, Juan J. Garcia Vallejo. Manuscript in preparation.

[^2]: Mage PL, Konecny AJ, Mair F. Measurement and prediction of unmixing-dependent spreading in spectral flow cytometry panels. bioRxiv [Preprint]. 2025. doi: 10.1101/2025.04.17.649396.
