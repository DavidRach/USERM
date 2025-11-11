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
#get Residual Obj
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
# fluors_selected = c(Sig_info$id[c(32:38,41:42,46:48,50:54,8,9,12,63)])#,58:62,8,9,12,14,18,24,25,26,63
fluors_selected = c(Sig_info$id[c(32:38,41:42,44,46:48,50:54,58,59,60,61,62,8,9,12,14,18,24,25,26,63)])#,58:62,8,9,12,14,18,24,25,26,63
Sig_info$id

fluors_selected = c(Sig_info$id[c(149,201,150,151,
                                  152,153,154,155,
                                  156,157,158,159,
                                  160,161,162,163,
                                  189,200,196,199,138)])#,58:62,8,9,12,14,18,24,25,26,63
# fluors_selected = c(Sig_info$id[c(149,201,150,151)])#,58:62,8,9,12,14,18,24,25,26,63
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
UsermObj$Scale_df$min = -2000
UsermObj$Scale_df$max = 5000
UsermObj$Scale_df$scale = "Linear"
UsermObj$Scale_df$cofactor = 10
PredOneSpread(Userm = UsermObj,population_id = c("V1"))
# PredOneSpread_update(Userm = UsermObj,population_id = c("V1"))

#create 3 populations with all zero fluroescence intensities
UsermObj$Intensity_mtx[,1] = 0
UsermObj$Intensity_mtx[,2] = UsermObj$Intensity_mtx[,1]
UsermObj$Intensity_mtx[,3] = UsermObj$Intensity_mtx[,1]
# assign 100, 500, 1000 to the PECy7 intensities of these 3 populations
UsermObj$Intensity_mtx[20,1] = 100
UsermObj$Intensity_mtx[20,2] = 500
UsermObj$Intensity_mtx[20,3] = 1000


UsermObj$Intensity_mtx[,2] = 200*c(1:length(UsermObj$fluors))
UsermObj$Intensity_mtx[,3] = 300*c(1:length(UsermObj$fluors))
names(UsermObj$Intensity_mtx) = c("V1","V2","V3")
PredMultipleSpread(Userm = UsermObj,population_ids = c("V1","V2","V3"))


# Step 4 matrics estimation and visualization
Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
# pdf(file = "E:/ResidualModel/CoefMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 1,min = 0,legend_name = "Coef",
        title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
# dev.off()

ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:31],
                  A = UsermObj$A,quiet = T)
pdf(file = "E:/ResidualModel/ssmMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = ssm,mincolor = "#B02B38",midcolor = "white", maxcolor = "#95ABDB",
        max = 2,mid = 0,min = -4,legend_name = "ss",
        title = "Spillover Spreading Matrix")
dev.off()
ss_output = check_ss(f_pos = "SCC_Cell_TCRVa24Ja18_BV785",
                     f_neg = "SCC_Bead_CD69_BV750",
                     SSM_fluor = UsermObj$fluors[1:31],
                     A = UsermObj$A,
                     custom_ssm_dir = "E:/MyFolder/ssm")
vis_ss_scatter(ss_output)


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

Spr1 = EstimateSpread(Userm = UsermObj,population_id = c("V1"))
Spr2 = EstimateSpread(Userm = UsermObj,population_id = c("V2"))
SpreadDistance_mtx = EstimateDistance(Spr1, Spr2)
Vis_Mtx(mat = SpreadDistance_mtx,mincolor = "darkred",midcolor = "white", maxcolor = "white",max = 1.5,mid = 1,min = 0,title = "Spread Distance matrix")


A = UsermObj$A
# pdf(file = "E:/ResidualModel/Signature matrix.pdf",width = 6,height = 12)
Vis_Mtx(mat = A,mincolor = "white",midcolor = "#D03E4C", maxcolor = "#B02B38",
        max = 1,mid = 0.5,min = 0,legend_name = "Signal",
        title = "Signature matrix")
# dev.off()

library(MASS)
A_pinv = ginv(A)
colnames(A_pinv) = rownames(A)
rownames(A_pinv) = colnames(A)
Vis_Mtx(mat = A_pinv,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 1,mid = 0,min = -1,legend_name = "Value",
        title = "Pseudo-inverse matrix")

A_pinv_t = t(A_pinv)
pdf(file = "E:/ResidualModel/transpose pinv matrix.pdf",width = 6,height = 12)
Vis_Mtx(mat = A_pinv_t,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 1,mid = 0,min = -1,legend_name = "Value",
        title = "transpose of Pseudo-inverse matrix")
dev.off()

pdf(file = "E:/ResidualModel/slop matrix.pdf",width = 12,height = 12)
slop_mtx = UsermObj$Res$SCC3_Cell_PECy7_CD4$slopMtx
Vis_Mtx(mat = slop_mtx,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 1,mid = 0,min = -1,legend_name = "beta",
        title = "slop matrix")
dev.off()

weight_mtx = A_pinv["SCC3_Cell_BV510_CD4",]%o%A_pinv["SCC3_Cell_BV510_CD4",]
pdf(file = "E:/ResidualModel/weight matrix.pdf",width = 12,height = 12)
Vis_Mtx(mat = weight_mtx,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 1,mid = 0,min = -1,legend_name = "weight",
        title = "weight matrix")
dev.off()


weighted_slope_mtx = weight_mtx * slop_mtx
pdf(file = "E:/ResidualModel/weight slope matrix.pdf",width = 12,height = 12)
Vis_Mtx(mat = weighted_slope_mtx,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 10,mid = 0,min = -10,legend_name = "value",
        title = "weighted slop matrix")
dev.off()
sum(weighted_slope_mtx)


ResObj = UsermObj$Res$SCC3_Cell_PECy7_CD4
ResObj$detectors
checkRes_covScatter(Res = ResObj,
                    detector1 = "B13-A",
                    detector2 = "B14-A")

checkSig_linePlot(id = "SCC3_Cell_PECy7_CD4")

#ssm
# ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:20],
#                   A = UsermObj$A,
#                   custom_ssm_dir = "E:/MyFolder/ssm")
ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:20],
                  A = UsermObj$A,quiet = T)
# A2 = UsermObj$A[-1,]
for(i in 1:nrow(UsermObj$A)){
  print(i)
  ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:20],
                    A = UsermObj$A[-i,],quiet = T)
  print(sum(abs(ssm[,8])))
}
Sig_info$id[155]
all_fluor = Sig_info$id[138:201]
fluors_selected = c(Sig_info$id[c(138,149,201,150,151,
                                  152,153,154,
                                  156,157,158,159,
                                  160,161,162,163,
                                  189,200,196,199,155)])
candidate_fluor = all_fluor[!(all_fluor %in% fluors_selected)]

for (new_fluor in candidate_fluor) {
  print(new_fluor)
  fluors_selected[21] = "SCC3_Cell_NFB530_CD4"#"SCC3_Cell_BUV805_CD4"#new_fluor

  Sig_mtx  = getSigMtx(ids = fluors_selected)
  UsermObj = CreateUserm(A = Sig_mtx)
  #add ResObj into UsermObj
  for (save_suf in colnames(Sig_mtx)) {
    ResObj = getRes(id = save_suf)
    UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
  }
  ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[2:21],
                    A = UsermObj$A,quiet = T)
  # print(sum(abs(ssm[,20])))
  print(max(ssm))
}

Sig_mtx  = getSigMtx(ids = fluors_selected)
UsermObj = CreateUserm(A = Sig_mtx)
#add ResObj into UsermObj
for (save_suf in colnames(Sig_mtx)) {
  ResObj = getRes(id = save_suf)
  UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
}

set.seed(123)

rand_sig_best = c()
A = UsermObj$A
ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[2:21],
                  A = A,quiet = T)
print(sum(abs(ssm[,20])))

A = cbind(A,runif(n=nrow(A), min=0, max=1))
colnames(A)[ncol(A)] = "pseudo"
ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[2:21],
                  A = A,quiet = T)
index_best = sum(abs(ssm[,20]))
ssm_best = ssm
print(index_best)
re_use = F
for (i in 1:1000) {
  if(!re_use){
    rand_sig = runif(n=nrow(A), min=0, max=0.1)
  }
  A_tmp = A
  A_tmp[,"pseudo"] = A_tmp[,"pseudo"] + rand_sig
  A_tmp[,"pseudo"] = A_tmp[,"pseudo"] / max(A_tmp[,"pseudo"])

  ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[2:21],
                    A = A_tmp,quiet = T)
  index = sum(abs(ssm[,20]))
  print(index)

  if (index < index_best){
    A = A_tmp
    index_best = index
    ssm_best = ssm
    re_use = T
  }else{
    re_use = F
  }

}

A[,"pseudo"]
Vis_Mtx(mat = ssm_best,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 10,mid = 0,min = -1,legend_name = "ss",
        title = "Spillover Spreading Matrix")


Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
# pdf(file = "E:/ResidualModel/CoefMtx.pdf",width = 10,height = 10)
Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
        max = 20,mid = 10,min = 0,legend_name = "Coef",
        title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
# dev.off()



Vis_Mtx(mat = ssm,mincolor = "#95ABDB",midcolor = "white", maxcolor = "#B02B38",
        max = 10,mid = 0,min = -1,legend_name = "ss",
        title = "Spillover Spreading Matrix")
max(ssm)
#check scatter for ss
ss_output = check_ss(f_pos = "SCC_Cell_CD16_BUV805",
                     f_neg = "SCC_Cell_CD19_BUV496",
                     SSM_fluor = UsermObj$fluors[1:20],
                     A = UsermObj$A,
                     custom_ssm_dir = "E:/MyFolder/ssm")
vis_ss_scatter(ss_output)
