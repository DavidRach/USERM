# devtools::install_github("xiangmingcai/USERM")
library(USERM)

# Step 1 querySig
Sig_info = querySig()
head(Sig_info)
fluor_to_check = Sig_info$id[3]
print(fluor_to_check)
checkSig_linePlot(id = fluor_to_check)

#get Residual Obj
#See \code{\link{CreateRes}} and \code{\link{SlopEstimation}} for detailed introduction of the structure of ResObj.

ResObj = getRes(id = fluor_to_check)
checkRes_slopMtx(Res = ResObj)
checkRes_interceptMtx(Res = ResObj)
print(ResObj$bin_mids)
checkRes_covMtx(Res = ResObj,bin=3)
ResObj$detectors
checkRes_covScatter(Res = ResObj,
                    detector1 = ResObj$detectors[1],
                    detector2 = ResObj$detectors[2])

# Step 2 Select fluors and create UsermObj
fluors_selected = c(Sig_info$id[c(32:38,41:42,46:48,50:54,8,9,12,63)])#,58:62,8,9,12,14,18,24,25,26,63
print(fluors_selected)
Sig_mtx  = getSigMtx(ids = fluors_selected)
dim(Sig_mtx)

UsermObj = CreateUserm(A = Sig_mtx)
#add ResObj into UsermObj
for (save_suf in colnames(Sig_mtx)) {
  ResObj = getRes(id = save_suf)
  UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
}

# Step 3 Make prediction (optional)
PredOneSpread(Userm = UsermObj,population_id = c("V1"))

UsermObj$Intensity_mtx[,2] = 200*c(1:length(UsermObj$fluors))
UsermObj$Intensity_mtx[,3] = 300*c(1:length(UsermObj$fluors))
names(UsermObj$Intensity_mtx) = c("P1","P2","P3")
PredMultipleSpread(Userm = UsermObj,population_ids = c("P1","P2","P3"))


# Step 4 matrics estimation and visualization
Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
# pdf(file = "E:/ResidualModel/CoefMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Coef",
        title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
# dev.off()

Hotspot_mtx = EstimateHotspotMtx(A = UsermObj$A)
# pdf(file = "E:/ResidualModel/HotspotMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = Hotspot_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Hotspot",
        title = "Hotspot matrix")
# dev.off()

Similarity_mtx = EstimateSimilarityMtx(A = UsermObj$A)
# pdf(file = "E:/ResidualModel/SimilarityMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = Similarity_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 1,mid = 0.8,min = 0,legend_name = "Cosine",
        title = "Cosine similarity matrix")
# dev.off()

Spr1 = EstimateSpread(Userm = UsermObj,population_id = c("P1"))
Spr2 = EstimateSpread(Userm = UsermObj,population_id = c("P2"))
SpreadDistance_mtx = EstimateDistance(Spr1, Spr2)
Vis_Mtx(mat = SpreadDistance_mtx,mincolor = "darkred",midcolor = "white", maxcolor = "white",max = 1.5,mid = 1,min = 0,title = "Spread Distance matrix")
