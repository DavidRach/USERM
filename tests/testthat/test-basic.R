test_that("querySig works", {
  # Step 1 querySig
  Sig_info = querySig()
  expect_equal(class(Sig_info), "data.frame")
  fluor_to_check = Sig_info$id[3]
  expect_equal(class(fluor_to_check), "character")
  expect_equal(class(checkSig_linePlot(id = fluor_to_check))[1], "gg")

  ResObj = getRes(id = fluor_to_check)
  expect_equal(class(ResObj), "list")
  expect_equal(class(checkRes_slopMtx(Res = ResObj))[1], "Heatmap")
  expect_equal(class(checkRes_interceptMtx(Res = ResObj))[1], "Heatmap")
  expect_equal(class(checkRes_covMtx(Res = ResObj,bin=3))[1], "Heatmap")
  expect_equal(class(checkRes_covScatter(Res = ResObj,
                                         detector1 = ResObj$detectors[1],
                                         detector2 = ResObj$detectors[2]))[1], "gg")
})

test_that("CreateUserm workflow", {
  # Step 1 querySig
  Sig_info = querySig()
  # Step 2 Select fluors and create UsermObj
  fluors_selected = c(Sig_info$id[c(32:34,63)])
  Sig_mtx  = getSigMtx(ids = fluors_selected)
  UsermObj = CreateUserm(A = Sig_mtx)
  expect_true(length(fluors_selected) > 0)
  expect_true(is.matrix(Sig_mtx))

  #add ResObj into UsermObj
  for (save_suf in colnames(Sig_mtx)) {
    ResObj = getRes(id = save_suf)
    UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
  }
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_LD_LDNIR876",new_name = "LiveDead")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD2_FITC",new_name = "CD2_FITC")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD3_BV510",new_name = "CD3_BV510")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_AF_AF",new_name = "AF")

  expect_type(UsermObj, "list")
  expect_type(ResObj, "list")
})

test_that("prediction works", {
  # Step 1 querySig
  Sig_info = querySig()
  # Step 2 Select fluors and create UsermObj
  fluors_selected = c(Sig_info$id[c(32:34,63)])
  Sig_mtx  = getSigMtx(ids = fluors_selected)
  UsermObj = CreateUserm(A = Sig_mtx)
  #add ResObj into UsermObj
  for (save_suf in colnames(Sig_mtx)) {
    ResObj = getRes(id = save_suf)
    UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
  }
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_LD_LDNIR876",new_name = "LiveDead")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD2_FITC",new_name = "CD2_FITC")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD3_BV510",new_name = "CD3_BV510")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_AF_AF",new_name = "AF")

  # Step 3 Make prediction (optional)
  UsermObj$Intensity_mtx[,2] = 200*c(1:length(UsermObj$fluors))
  UsermObj$Intensity_mtx[,3] = 300*c(1:length(UsermObj$fluors))
  names(UsermObj$Intensity_mtx) = c("P1","P2","P3")

  expect_type(PredOneSpread(Userm = UsermObj,population_id = c("P1")), "list")
  expect_type(PredMultipleSpread(Userm = UsermObj,population_ids = c("P1","P2","P3")), "list")



})

test_that("visualization workflow", {
  # Step 1 querySig
  Sig_info = querySig()
  # Step 2 Select fluors and create UsermObj
  fluors_selected = c(Sig_info$id[c(32:34,63)])
  Sig_mtx  = getSigMtx(ids = fluors_selected)
  UsermObj = CreateUserm(A = Sig_mtx)
  expect_true(length(fluors_selected) > 0)
  expect_true(is.matrix(Sig_mtx))

  #add ResObj into UsermObj
  for (save_suf in colnames(Sig_mtx)) {
    ResObj = getRes(id = save_suf)
    UsermObj = AddRes2Userm(Res = ResObj, Userm = UsermObj)
  }
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_LD_LDNIR876",new_name = "LiveDead")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD2_FITC",new_name = "CD2_FITC")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_CD3_BV510",new_name = "CD3_BV510")
  UsermObj = RenameFluor(Userm = UsermObj,raw_name = "SCC_Cell_AF_AF",new_name = "AF")
  # Step 3 Make prediction (optional)
  UsermObj$Intensity_mtx[,2] = 200*c(1:length(UsermObj$fluors))
  UsermObj$Intensity_mtx[,3] = 300*c(1:length(UsermObj$fluors))
  names(UsermObj$Intensity_mtx) = c("P1","P2","P3")


  # Step 4 matrics estimation and visualization
  Coef_mtx = EstimateCoefMtx(Userm = UsermObj)
  expect_equal(is.matrix(Coef_mtx), TRUE)
  Vis_Mtx(mat = Coef_mtx,mincolor = "white",midcolor = "white", maxcolor = "#95ABDB",
          max = 2,mid = 1,min = 0,legend_name = "Coef",
          title = "Coefficient of residual model matrix (Coef Matrix) (row spread into column)")
  ssm = EstimateSSM(SSM_fluor = UsermObj$fluors[1:3],
                    Userm = UsermObj,
                    A = UsermObj$A,quiet = T)
  expect_equal(is.matrix(ssm), TRUE)
  ss_output = check_ss(f_pos = "CD2_FITC",
                       f_neg = "CD3_BV510",
                       SSM_fluor = UsermObj$fluors[1:3],
                       A = UsermObj$A,
                       Userm = UsermObj)
  expect_equal(class(ss_output), "list")
  expect_equal(class(vis_ss_scatter(ss_output))[1], "gg")
  Hotspot_mtx = EstimateHotspotMtx(A = UsermObj$A)
  expect_equal(is.matrix(Hotspot_mtx), TRUE)
  Similarity_mtx = EstimateSimilarityMtx(A = UsermObj$A)
  expect_equal(is.matrix(Similarity_mtx), TRUE)

  Spr1 = EstimateSpread(Userm = UsermObj,population_id = c("P1"))
  Spr2 = EstimateSpread(Userm = UsermObj,population_id = c("P2"))
  SpreadDistance_mtx = EstimateDistance(Spr1, Spr2)
  expect_equal(class(Spr1), "list")
  expect_equal(is.matrix(SpreadDistance_mtx), TRUE)
})
