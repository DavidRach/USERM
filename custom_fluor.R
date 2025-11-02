devtools::install_github("xiangmingcai/USERM@dev")
# devtools::install_github("xiangmingcai/USERM")
library(USERM)
library(flowCore)
# devtools::install_github("xiangmingcai/GateData")
library(GateData)
library(dplyr)
library(ggplot2)

custom_dir = "E:/MyFolder"
dir.create(paste0(custom_dir,"/sig"))
dir.create(paste0(custom_dir,"/res"))

#### Part 1 Prepare signature ####
#step 1.1 read scc fcs
data = read.FCS("E:/Data/SCC_Cell_CD2_SB780 Run 1 20251024102726.fcs")
# data = read.FCS("E:/Data/unstained_1 Run 1 20251024100712.fcs")

print(head(data@parameters@data$desc))
print(head(data@parameters@data$name))
desc = data@parameters@data$desc # "desc" or "name" are used for different insruments. Ues the one with correct detector names.

data = exprs(data)
data = as.data.frame(data)
colnames(data) = desc
head(data)

# step 1.2 set scc fcs info
peak_channel = '405nm - 770/LP-A' #use public resource (e.g. fluorofinder) to find the peak channel of your SCC.
save_suf = "SCCcustom_Cell_CD2_SB780"
# save_suf = "SCCcustom_Cell_AF_AF"
#it is recommended to follow the naming pattern.
#There are 4 elements in the name, which are seperated by underscore.
#The first is "batch" element. The preprocessed scc use "SCC1", "SCC2", and so on to annotate SCCs acquired in distinct batches.
#You can use "SCCcustom" to distinguish custom fluorescence and preprocessed fluorescence. You can also use others that make sense to you.
#The second is "type" element. Now We only have "Cell" and "Bead".
#The third is "PrimaryName". Normally we put target of the fluorescence antibody here.
#The last is "SecondaryName", Normally we put fluroescence of the antibody here.
#Please avoid using space or underscore in these elements.
PrimaryName = "CD2"
SecondaryName = "SB780"
# PrimaryName = "AF"
# SecondaryName = "AF"

#step 1.3 gate positive popualtion and negative population
#In this step, we need to find positive and negative populations. You can use the GateData R package or other packages that works for you. Here we show how to gate a subset population with the GateData
#For detailed instruction of GateData, please refer to https://github.com/xiangmingcai/GateData

# 1. sample cells.
#You may modify code here to sample specific set for rare markers.
data$gate0 = TRUE
data = sample_n(data, 20000)

# 2. gate singlets if need
#You can change the FSC and SSC names for specific instrument.
# use colnames(data) to find available paramters
# colnames(data)
gate1<-PolygonGating(df=data, x_col= "488 FSC-A", y_col= "488 FSC-H", feature_col= "488 FSC-A",
                     parentgate_col= "gate0", newgate_col= "gate1",canvas_width=800, canvas_height=400,
                     title_text = "Gate singlets")
data <-GateDecider(gate = gate1, df = data)
gate2<-PolygonGating(df=data, x_col= "488 SSC-A", y_col= "488 SSC-H", feature_col= "488 SSC-A",
                     parentgate_col= "gate1", newgate_col= "gate2",canvas_width=800, canvas_height=400,
                     title_text = "Gate singlets")
data <-GateDecider(gate = gate2, df = data)

# 3. gate target population (e.g. lymphocytes)
gate3<-PolygonGating(df=data, x_col= "488 SSC-A", y_col= "488 FSC-A", feature_col= "488 SSC-A",
                     parentgate_col= "gate2", newgate_col= "gate3",canvas_width=1500, canvas_height=1500,
                     title_text = "Gate target population")
data <-GateDecider(gate = gate3, df = data)

# 4. gate positive and negative populations.
#For autofluorescence, randomly gate a negative populaiton and assign zero to it.
data = data[data$gate3,]
gate_pos<-PolygonGating(df=data, x_col= peak_channel, y_col= "488 SSC-A", feature_col= "488 SSC-A",
                        parentgate_col= "gate3", newgate_col= "gate_pos",canvas_width=800, canvas_height=400,
                        title_text = "Gate positive population")
data <-GateDecider(gate = gate_pos, df = data)
gate_neg<-PolygonGating(df=data, x_col= peak_channel, y_col= "488 SSC-A", feature_col= "488 SSC-A",
                        parentgate_col= "gate3", newgate_col= "gate_neg",canvas_width=800, canvas_height=400,
                        title_text = "Gate negative population")
data <-GateDecider(gate = gate_neg, df = data)

table(data$gate_pos)
table(data$gate_neg)
data_pos = data[data$gate_pos,]
data_neg = data[data$gate_neg,]

# For autofluorescence, assign zero to data_neg
# data_neg[,] <- 0

# #you may save the gates if needed.
# saveRDS(list(gate1,gate2,gate3,gate_pos,gate_neg), file = paste0(tmp_dir,"/",save_suf,"_gatelist.rds"))
# saveRDS(c("gate1","gate2","gate3","gate_pos","gate_neg"), file = paste0(tmp_dir,"/",save_suf,"_gatenames.rds"))

#step 1.4 extract signatures
#Select detectors. Please select only the detectors for unmixing
colnames(data)
cols = colnames(data)[seq(from=28, to=128, by=2)]
print(cols)

Sig = ExtractSig(df_pos = data_pos,df_neg = data_neg,
                 cols = cols, method = "median",
                 PrimaryName = PrimaryName,
                 SecondaryName = SecondaryName,
                 id = save_suf,
                 instrument = "Xenith",
                 Source = "YourName",
                 Note = NA)

#you may briefly check the signature
df <- data.frame(Detector = names(Sig$Signature),Signal = as.numeric(Sig$Signature))
ggplot(df, aes(x = Detector, y = Signal, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = save_suf,
       x = "Detector",
       y = "Signal")

#step 1.5 add signatures to Custom_Sig_list
# 1. Read in Custom_Sig_list
#For first SCC (that means you do not have a Custom_Sig_list yet)
Custom_Sig_list = list()
#If you already have a Custom_Sig_list
Custom_Sig_list = readRDS(paste0(custom_dir,"/sig/Custom_Sig_list.rds"))

# 2. add signatures to Custom_Sig_list
Custom_Sig_list[[length(Custom_Sig_list)+1]] = Sig
names(Custom_Sig_list)[length(Custom_Sig_list)] = Custom_Sig_list[[length(Custom_Sig_list)]]$id

# You may delete a certain Sig from the Custom_Sig_list using this code:
# Custom_Sig_list[[length(Custom_Sig_list)]] = NULL

# 3. save updated Custom_Sig_list
saveRDS(Custom_Sig_list, file = paste0(custom_dir,"/sig/Custom_Sig_list.rds"))



#### Part 2 Prepare res obj ####
#Now you should have at least one SCC signature and one Autofluroescence signature.

#step 2.1 read Custom_Sig_list and prepare Sig_mtx
Custom_Sig_list = readRDS(paste0(custom_dir,"/sig/Custom_Sig_list.rds"))

Sig_info = querySig(Sig_list = Custom_Sig_list)
fluors_selected = c(Sig_info$id[c(1,2)])
Sig_mtx  = getSigMtx(ids = fluors_selected, Sig_list = Custom_Sig_list)

#step 2.2 read scc fcs
data = read.FCS("E:/Data/SCC_Cell_CD2_SB780 Run 1 20251024102726.fcs")
# data = read.FCS("E:/Data/unstained_1 Run 1 20251024100712.fcs")

print(head(data@parameters@data$desc))
print(head(data@parameters@data$name))
desc = data@parameters@data$desc # "desc" or "name" are used for different insruments. Ues the one with correct detector names.

data = exprs(data)
data = as.data.frame(data)
colnames(data) = desc
head(data)

# step 2.3 set target and AF info
save_suf = "SCCcustom_Cell_CD2_SB780"
idx_target_fluor = 1
idx_AF_fluor = 2

# For AF
# save_suf = "SCCcustom_Cell_AF_AF"
# idx_target_fluor = 2
# idx_AF_fluor = 1#For AF, use any other scc here

#check if the idx are correct
colnames(Sig_mtx)[idx_target_fluor]
colnames(Sig_mtx)[idx_AF_fluor]

# step 2.4 gate signlets and target population
#Here we want to use cells that falls into gate as many as possible. So, we use a sampled subset to make gates and subset from all cells.

# 1. sample cells.
#You may modify code here to sample specific set for rare markers.
data$gate0 = TRUE
data_backup = data

data = sample_n(data, 20000)

# 2. gate singlets if need
#You can change the FSC and SSC names for specific instrument.
# use colnames(data) to find available paramters
# colnames(data)
gate1<-PolygonGating(df=data, x_col= "488 FSC-A", y_col= "488 FSC-H", feature_col= "488 FSC-A",
                     parentgate_col= "gate0", newgate_col= "gate1",canvas_width=800, canvas_height=400,
                     title_text = "Gate singlets")
data <-GateDecider(gate = gate1, df = data)
gate2<-PolygonGating(df=data, x_col= "488 SSC-A", y_col= "488 SSC-H", feature_col= "488 SSC-A",
                     parentgate_col= "gate1", newgate_col= "gate2",canvas_width=800, canvas_height=400,
                     title_text = "Gate singlets")
data <-GateDecider(gate = gate2, df = data)

# 3. gate target population (e.g. lymphocytes)
gate3<-PolygonGating(df=data, x_col= "488 SSC-A", y_col= "488 FSC-A", feature_col= "488 SSC-A",
                     parentgate_col= "gate2", newgate_col= "gate3",canvas_width=1500, canvas_height=1500,
                     title_text = "Gate target population")
data <-GateDecider(gate = gate3, df = data)

# 4. subset with all cells
data = data_backup
data <-GateDecider(gate = gate1, df = data)
data <-GateDecider(gate = gate2, df = data)
data <-GateDecider(gate = gate3, df = data)
data = data[data$gate3,]

#step 2.5 unmix with target_fluor and AF
data = data[,rownames(Sig_mtx)] #cells x detectors
R = t(as.matrix(data)) # detectors x cells
A = Sig_mtx[,c(idx_target_fluor,idx_AF_fluor), drop = FALSE] # detectors x fluors
A_pinv = ginv(A)
B =  A_pinv %*% R
t_B = t(B)
colnames(t_B) = c(colnames(Sig_mtx)[idx_target_fluor],colnames(Sig_mtx)[idx_AF_fluor])
data = cbind(data,t_B)

# 2.6 gate normal samples

# You may remove some outliers for easier gating
# data = data[(data[,colnames(Sig_mtx)[idx_target_fluor]] <
#                quantile(data[,colnames(Sig_mtx)[idx_target_fluor]],0.999)),]
# data = data[(data[,colnames(Sig_mtx)[idx_target_fluor]] >
#                quantile(data[,colnames(Sig_mtx)[idx_target_fluor]],0.01)),]

data_backup = data
data = sample_n(data, 20000)
gate_normal<-PolygonGating(df=data, x_col= colnames(Sig_mtx)[idx_target_fluor],
                           y_col= colnames(Sig_mtx)[idx_AF_fluor],
                           feature_col= colnames(Sig_mtx)[idx_target_fluor],
                           parentgate_col= "gate0", newgate_col= "gate_normal",canvas_width=800, canvas_height=400)
data = data_backup
data = GateDecider(gate = gate_normal, df = data)
data = data[data$gate_normal,rownames(Sig_mtx)]

#you may save gate_normal is you want
# saveRDS(list(gate_normal), file = paste0(tmp_dir,"/",save_suf,"_gatelist.rds"))
# saveRDS(c("gate_normal"), file = paste0(tmp_dir,"/",save_suf,"_gatenames.rds"))

# 2.7 calculate and save residual model
R = as.matrix(data)
A_Target = Sig_mtx[,c(idx_target_fluor), drop = FALSE]
A_AF = Sig_mtx[,c(idx_AF_fluor), drop = FALSE]

ResObj = CreateRes(id = colnames(Sig_mtx)[idx_target_fluor], R = R, A_Target = A_Target, A_AF = A_AF)

ResObj = SlopEstimation(Res = ResObj, bin_num = 30)

#You may briefly visualize the slope matrix
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(ResObj$slopMtx,col = col_fun, cluster_rows = FALSE,cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", ResObj$slopMtx[i, j]), x, y, gp = gpar(fontsize = 6))
        },column_title = colnames(Sig_mtx)[idx_target_fluor],name = "beta")

saveRDS(ResObj, file = paste0(custom_dir,"/res","/ResObj_",save_suf,".rds"))





#### Part 3 Use custom res obj and preprocessed res obj together ####
# tbd



