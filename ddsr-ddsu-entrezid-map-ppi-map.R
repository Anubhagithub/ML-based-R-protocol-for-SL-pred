#map the gene names to entrez ids
start = proc.time()
library(doParallel)
library(foreach)
n.cores <- parallel::detectCores() - 5
c1 <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)
doParallel::registerDoParallel(c1)
# ddsr = read.csv("DD-SR.csv")
# dusr = read.csv("DU-SR.csv")
# synlegg = read.csv("synlegg_alltissue_ppi.csv")
# synlegg2 = synlegg[,c(1,2)]
# synlegg2 = unique(synlegg2)
# entrez = read.csv("unid_genename_entrezid_withoutspace.csv")
# entrez2 = entrez[,c(2,4)]
# entrez3 = na.omit(entrez2)
# entrez_cp = entrez3
# sv_cp = synlegg2
# result = c()
# result = c()
# result = foreach(x = 1:nrow(sv_cp),.packages='doParallel',.combine = 'rbind') %dopar% {
#   foreach(y = 1:nrow(entrez_cp), .packages='doParallel',.combine = 'rbind') %do%
#     {
#       if(grepl(sv_cp[x,1], entrez_cp[y,1]) && nchar(sv_cp[x,1]) == nchar(entrez_cp[y,1]))
#       {
#         foreach(z = 1:nrow(entrez_cp), .packages='doParallel',.combine = 'rbind') %do%
#           {
#             if(grepl(sv_cp[x,2], entrez_cp[z,1]) && nchar(sv_cp[x,2]) == nchar(entrez_cp[z,1]))
#             {
#               c(entrez_cp[y,2], entrez_cp[z,2])
#               #list2 = rbind(list2, c(net_cp[y,], net_cp[z,]))
#               #c(paste(c(net_cp[y,], net_cp[z,])))
#               #break
#             }
#           }
#       }
#     }
# }
# final = as.data.frame(result)
# end = proc.time()
# time_taken = end - start
# stopCluster(c1)
# write.csv(result, "synlegg_entrezid.csv", row.names = F, quote = F)
#calculate the nodewise netprop
#copy the above file generated to nikola
#map svs to ppi
# dusr_entrez = read.csv("dusr_entrezid.csv")
# dusr_entrez2 = na.omit(dusr_entrez)
synlegg_entrez = read.csv("synlegg_entrezid.csv")
anyNA(synlegg_entrez)
net = read.csv("nodewise_netprop_ppi.csv")
not_cp = synlegg_entrez
net_cp = net
result = c()
result = foreach(x = 1:nrow(not_cp),.packages='doParallel',.combine = 'rbind') %dopar% {
  foreach(y = 1:nrow(net_cp), .packages='doParallel',.combine = 'rbind') %do%
    {
      if(grepl(not_cp[x,1], net_cp[y,1]) && nchar(not_cp[x,1]) == nchar(net_cp[y,1]))
      {
        foreach(z = 1:nrow(net_cp), .packages='doParallel',.combine = 'rbind') %do%
          {
            if(grepl(not_cp[x,2], net_cp[z,1]) && nchar(not_cp[x,2]) == nchar(net_cp[z,1]))
            {
              c(net_cp[y,], net_cp[z,])
              #list2 = rbind(list2, c(net_cp[y,], net_cp[z,]))
              #c(paste(c(net_cp[y,], net_cp[z,])))
              #break
            }
          }
      }
    }
}
final = as.data.frame(result)
end = proc.time()
time_taken = end - start
stopCluster(c1)
write.csv(result, "synlegg_nodewise_netprop.csv", quote = F, row.names = F)
# # #calculate pairwise properties
# # inp2 = read.csv("ddsr_nodewise_netprop.csv")
# inp2 = read.csv("dusr_nodewise_netprop.csv")
inp2 = read.csv("synlegg_nodewise_netprop.csv")
inp3 = inp2[, c(1,16,2,17,3,18,4,19,5,20,6,21,7,22,8,23,9,24,10,25,11,26,12,27,13,28,14,29,15,30)]
inp3$names = as.integer(inp3$names)
inp3$names.1 = as.integer(inp3$names.1)
inp3$degree = as.integer(inp3$degree)
inp3$degree.1 = as.integer(inp3$degree.1)
inp3$betwn = as.integer(inp3$betwn)
inp3$betwn.1 = as.integer(inp3$betwn.1)
inp3$close = as.integer(inp3$close)
inp3$close.1 = as.integer(inp3$close.1)
inp3$core = as.integer(inp3$core)
inp3$core.1 = as.integer(inp3$core.1)
inp3$constraint = as.integer(inp3$constraint)
inp3$constraint.1 = as.integer(inp3$constraint.1)
inp3$ecc = as.integer(inp3$ecc)
inp3$ecc.1 = as.integer(inp3$ecc.1)
inp3$eigen_cen = as.integer(inp3$eigen_cen)
inp3$eigen_cen.1 = as.integer(inp3$eigen_cen.1)
inp3$hub_score = as.integer(inp3$hub_score)
inp3$hub_score.1 = as.integer(inp3$hub_score.1)
inp3$neighbor1 = as.integer(inp3$neighbor1)
inp3$neighbor1.1 = as.integer(inp3$neighbor1.1)
inp3$neighbor2 = as.integer(inp3$neighbor2)
inp3$neighbor2.1 = as.integer(inp3$neighbor2.1)
inp3$neighbor3 = as.integer(inp3$neighbor3)
inp3$neighbor3.1 = as.integer(inp3$neighbor3.1)
inp3$neighbor4 = as.integer(inp3$neighbor4)
inp3$neighbor4.1 = as.integer(inp3$neighbor4.1)
inp3$neighbor5 = as.integer(inp3$neighbor5)
inp3$neighbor5.1 = as.integer(inp3$neighbor5.1)
inp3$neighbor6 = as.integer(inp3$neighbor6)
inp3$neighbor6.1 = as.integer(inp3$neighbor6.1)
#calculating avg for all properties
inp3$avgdeg = (inp3$degree+inp3$degree.1)/2
inp3$avgbet = (inp3$betwn+inp3$betwn.1)/2
inp3$avgclose = (inp3$close+inp3$close.1)/2
inp3$avgcore = (inp3$core+inp3$core.1)/2
inp3$avgcons = (inp3$constraint+inp3$constraint.1)/2
inp3$avgecc = (inp3$ecc+inp3$ecc.1)/2
inp3$avgeigen_cen = (inp3$eigen_cen+inp3$eigen_cen.1)/2
inp3$avghub_score = (inp3$hub_score+inp3$hub_score.1)/2
inp3$avgneighbor1 = (inp3$neighbor1+inp3$neighbor1.1)/2
inp3$avgneighbor2 = (inp3$neighbor2+inp3$neighbor2.1)/2
inp3$avgneighbor3 = (inp3$neighbor3+inp3$neighbor3.1)/2
inp3$avgneighbor4 = (inp3$neighbor4+inp3$neighbor4.1)/2
inp3$avgneighbor5 = (inp3$neighbor5+inp3$neighbor5.1)/2
inp3$avgneighbor6 = (inp3$neighbor6+inp3$neighbor6.1)/2
inp3cp = inp3[, c(1:2,31:44)]
inp4 = na.omit(inp3cp)
dim(inp4)
colnames(inp4)[1] = "gene1"
colnames(inp4)[2] = "gene2"
inp4$gene1 = as.factor(inp4$gene1)
inp4$gene2 = as.factor(inp4$gene2)
for (i in 1:length(inp4$gene1)){
  V = inp4$gene1[i]
  U = inp4$gene2[i]
  inp4$shortest_path[i] <- distances(Network2, v = V, to = U)
}
#inp4 = inp4[, -18]
inp4$cohesion = 0
for (i in 1:length(inp4$gene1)){
  if (inp4$shortest_path[i] > 1) {
    V = inp4$gene1[i]
    U = inp4$gene2[i]
    print(i)
    inp4$cohesion[i] <- vertex.connectivity(Network2, source = V,target = U)
  }
}
inp4[inp4$shortest_path<=1,]$cohesion=NA
#inp4 = inp4[, -19]
inp4 = na.omit(inp4)
inp4$adhesion = 0
for (i in 1:length(inp4$gene1)){
  V = inp4$gene1[i]
  U = inp4$gene2[i]
  print(i)
  inp4$adhesion[i] <- edge.connectivity(Network2, source = V, target = U)
}
inp4 = na.omit(inp4)
write.csv(inp4,"synlegg_pairwisenetprop.csv", quote = F, row.names = F)

#prepare the input for tseting the model
one = read.csv("synlegg_pairwisenetprop.csv")
one$gi = "SL"
two = read.csv("ddsr_pairwisenetprop.csv")
two$gi = "SV"
three = read.csv("dusr_pairwisenetprop.csv")
three$gi = "SV"
four = read.csv("/home/nikola/biogrid_sl_sv_not/idmap/entrezid/NOT/98_not_nodewise_net_prop.csv")
four$gi = "NOT"
valdt = rbind(one,two,three,four)
anyNA(valdt)
valdt = valdt[which(valdt$shortest_path != Inf),]
write.csv(valdt, "validation_data_synlegg_ddsr_dusr_not98.csv", row.names = F,quote = F)
#is there ant overlap between ddsr and dusr
a = two
b = three
a = two[,c(1,2)]
b = three[,c(1,2)]
a$pair = paste(a$gene1, a$gene2, sep = "_")
b$pair = paste(b$gene1, b$gene2, sep = "_")
acp = as.data.frame(a[,c(3)])
bcp = as.data.frame(b[,c(3)])
colnames(acp)[1] = "pair"
colnames(bcp)[1] = "pair"
intr_acp_bcp = intersect(acp, bcp)
#read the file generated
valdt = read.csv("/home/nikola/sldatasets/validation_data_synlegg_ddsr_dusr_not98.csv")
valdt = valdt[which(valdt$shortest_path != Inf),]
valdt2 = valdt[,c(3:20)]
valdt2$gi = as.factor(valdt2$gi)
table(valdt2$gi)
#training file
cgidata = read.csv("/home/nikola/biogrid_sl_sv_not/idmap/work_after_4aug_drc/cgidb_olddata_netprop.csv")
cgidata2 = cgidata[which(cgidata$shortest_path != Inf),]
cgidata2$gi = as.factor(cgidata2$gi)
table(cgidata2$gi)
biodata = read.csv("/home/nikola/sldatasets/biogrid_netprop.csv")
biodata2 = biodata[which(biodata$shortest_path != Inf),]
biodata2$gi = as.factor(biodata2$gi)
biodata2 = biodata2[,c(3:20)]
not1 = read.csv("/home/nikola/biogrid_sl_sv_not/idmap/entrezid/NOT/100_not_nodewise_net_prop.csv")
not1 = not1[which(not1$shortest_path != Inf),]
not1$gi = "NOT"
not1$gi = as.factor(not1$gi)
not1 = not1[,c(3:20)]
datatotrain = unique(rbind(cgidata2,biodata2))
slkasldb = read.csv("/home/nikola/sldatasets/sl.sldb_pairwisenetprop.csv")
anyNA(slkasldb)
slkasldb$gi = "SL"
slkasldb = slkasldb[,c(3:20)]
svkasvdr = read.csv("/home/nikola/sldatasets/svdr_pairwise_netprop.csv")
anyNA(svkasvdr)
svkasvdr = svkasvdr[which(svkasvdr$shortest_path != Inf),]
not2 = read.csv("/home/nikola/biogrid_sl_sv_not/idmap/entrezid/NOT/99_not_nodewise_net_prop.csv")
not2 = not2[which(not2$shortest_path != Inf),]
not2$gi = "NOT"
not2$gi = as.factor(not2$gi)
not2 = not2[,c(3:20)]
sldbsvdrnotnew = unique(rbind(slkasldb,svkasvdr,not2))
sldbsvdrnotnew = sldbsvdrnotnew[which(sldbsvdrnotnew$shortest_path != Inf),]
sldbsvdrnotnew$gi = as.factor(sldbsvdrnotnew$gi)
revise_train = unique(rbind(datatotrain,sldbsvdrnotnew))
extra = unique(union_all(datatotrain,sldbsvdrnotnew))
#newdata = rbind(slofsldb,sv.svdr2,not_ours2)
#newdata2 = newdata[which(newdata$shortest_path != Inf),]
#newdata2$gi = as.factor(newdata2$gi)
#check overlap of train and validation
train.intr = intersect(datatotrain, sldbsvdrnotnew)
train.val.intr = intersect(revise_train, valdt2)
diff = setdiff(revise_train, train.val.intr)
intr_crosscheck = intersect(diff, valdt2)
write.csv(diff, "cgi.bio.sldb.svdr.not.nooverlapwithvaldation.csv", row.names = F, quote = F)
#now assign train and test 
valdt.cp = read.csv("/home/nikola/brain-ppi-revise/valdt.csv")
valdt2 = valdt.cp[,c(4:21)]
valdt2$gi = as.factor(valdt2$gi)
diff.bal.cp = read.csv("/home/nikola/brain-ppi-revise/diff.bal.csv")
diff.bal = diff.bal.cp[,c(2:19)]
diff.bal$gi = as.factor(diff.bal$gi)
# main_fltr2 = diff
# index = sample(2, nrow(main_fltr2), replace = TRUE, prob = c(0.60, 0.40))
# training = main_fltr2[index == 1,]
# testing = main_fltr2[index == 2,]
set.seed(123)
training = diff.bal
testing = valdt2
library(randomForest)
model1 = randomForest(gi ~ . , data = training, importance = T)
prediction = predict(model1, newdata = testing[-18])
predict_gi = predict(model1, newdata = testing)
testing$predict_gi = predict_gi
cnf_mat = table(testing$gi, testing$predict_gi)
cnf_mat
accuracy = sum(diag(cnf_mat)/sum(cnf_mat))
accuracy2 = sum(diag(cnf_mat)/sum(cnf_mat))
training = diff
testing = valdt2
diff.bal = DMwR::SMOTE(gi ~ ., diff, perc.under = 200)
write.csv(diff.bal, "diff.bal.csv", row.names = T, quote = F)
varImpPlot(model1, sort = T, n.var = 17, main = "Discriminatory variables")
importance(model1)
varUsed(model1)
pred_score = predict(model1,testing, type = "vote", norm.vote = T)
head(pred_score)
istheresvintest = testing[which(testing$predict_gi == "SV"),]
#test scored
#https://rpubs.com/jkylearmstrong/RF_Imputation_Multi_class
probs <- predict(model1, testing, 'prob')
class <- predict(model1, testing)
TEST.scored <- cbind(testing, probs, class)
library('yardstick')
cm <- conf_mat(TEST.scored, truth = gi, class)
library('ggplot2')
ggplot(summary(cm), aes(x=.metric, y=.estimate)) + 
  geom_bar(stat="identity", fill = "#FF62BC") + 
  coord_flip()
#next step is feature selection
#we already know about the varImpPlot of Random forest now we use boruta to better
#understand the discriminatory features
diff.sub = diff[,c(2,18)]
valdt2.sub = valdt2[, c(2,18)]
training = diff.sub
testing = valdt2.sub
#error rate
head(model1$err.rate)
print(model1)
# Call:
#   randomForest(formula = gi ~ ., data = training) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 4
# 
# OOB estimate of  error rate: 13.85%
# Confusion matrix:
#   NOT    SL   SV class.error
# NOT 22473  2182    2  0.08857525
# SL   1737 12694 1020  0.17843505
# SV     64   933 1780  0.35902053
#trying to understand the error rate
library(ROCR)
OOB.votes <- predict(model1,x,type="prob")
OOB.pred <- OOB.votes[,2];

pred.obj <- prediction(OOB.pred,y)

RP.perf <- performance(pred.obj, "rec","prec");
plot (RP.perf);

ROC.perf <- performance(pred.obj, "fpr","tpr");
plot (ROC.perf);

plot  (RP.perf@alpha.values[[1]],RP.perf@x.values[[1]]);
lines (RP.perf@alpha.values[[1]],RP.perf@y.values[[1]]);
lines (ROC.perf@alpha.values[[1]],ROC.perf@x.values[[1]]);
#to find the correlation between the different features of the model
#p value and correlation
library(Hmisc)
library(corrplot)
corr.data = diff[,-c(18)]
cor_5 <- rcorr(as.matrix(corr.data))
M <- cor_5$r
p_mat <- cor_5$P
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p_mat, sig.level = 0.01)
#corr and p value for nodewise properties
cor_5 <- rcorr(as.matrix(netprop[, c(2:15)]))
anyNA(netprop)
M <- cor_5$r
p_mat <- cor_5$P
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p_mat, sig.level = 0.01)
sub_netprop = netprop[, c(1:18)]
sub_netprop = netprop[,c(2:11)]
cor_5 <- rcorr(as.matrix(sub_netprop[,]))
M <- cor_5$r
p_mat <- cor_5$P
corrplot(M, type = "upper", order = "hclust", 
         p.mat = p_mat, sig.level = 0.01, insig = "blank")
corrplot(M,  order = "hclust", 
         p.mat = p_mat, sig.level = 0.01, insig = "blank")
corrplot(M, method = 'number')
#boruta
diff.bal.cp = read.csv("/home/nikola/brain-ppi-revise/diff.bal.csv")
diff.bal = diff.bal.cp[,c(2:19)]
diff.bal$gi = as.factor(diff.bal$gi)
library(Boruta)
set.seed(111)
boruta.data = diff.bal
boruta.data$random<-sample(3, size = nrow(boruta.data), replace = TRUE)
boruta.data = boruta.data[,c(1:17,19,18)]
boruta <- Boruta(gi ~ ., data = boruta.data, doTrace = 2, maxRuns = 500)
print(boruta)
#actual plot
par(mar=c(5,4,0.28,1))
plot(boruta, add = T, las = 2, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"), 
     sort = TRUE, whichShadow = c(TRUE, TRUE, TRUE), col = NULL, 
     xlab = "Attributes", ylab = "Importance", cex.lab=1.2, col.lab = "red")
par(mar=c(6,4,0.15,1))
plot(boruta, las = 2, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"), 
     sort = TRUE, whichShadow = c(TRUE, TRUE, TRUE), col = NULL, 
     xlab = "Attributes", ylab = "Importance", cex.lab=1.6, col.lab = "red")
plot(boruta, las = 2, add = T, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"))
plot(boruta, las = 2, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"))
plotImpHistory(boruta)
library(ggplot2)
tiff("boruta.tiff", height = 40, width = 30, res = 100)
plot(boruta, add = T, las = 2, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"))
dev.off()
#trying to change boruta color
plot(boruta, add = TRUE, las = 2, cex.axis = 0.8,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"),
     xlab = substitute(paste(bold('X Label'))),
     ylab = substitute(paste(bold('Y Label'))))
plot(boruta, add = T, las = 2, cex.axis = 0.7,colCode = c("#00BC59", "yellow", "#FF6C90", "#00BBDB"), 
     sort = TRUE, whichShadow = c(TRUE, TRUE, TRUE), col = NULL, 
     xlab = "Attributes", ylab = "Importance", cex.lab=3.5)
#plot to show the how often the variables have been used for model1
Varused <- c(46386, 69953, 38227, 41092, 41982, 31757, 36987, 37451, 46199, 70972,
69615, 68860, 70228, 71286, 34851, 33019, 34332)
Varused = c(46273, 70521, 38494, 40941, 42183, 31818, 37395, 37114, 46318, 70746,
69415, 69556, 70452, 70862, 34778, 32866, 34162)
Network_prop <- c(rep("avgdeg" , 1) ,rep("avgbet" , 1) , rep("avgclose" , 1) , rep("avgcore" , 1), rep("avgcons", 1),
rep("avgecc" , 1), rep("avgeigcen" , 1), rep("avghub" , 1), rep("avgneighbor1" , 1), rep("avgneighbor2" , 1),
rep("avgneighbor3" , 1), rep("avgneighbor4" , 1), rep("avgneighbor5" , 1), rep("avgneighbor6" , 1), rep("shortest_path" , 1),
rep("cohesion" , 1), rep("adhesion" , 1))
data <- data.frame(Varused, Network_prop)    
Network_property = reorder(Network_prop,-Varused)
q = ggplot(data, las = 1, cex.axis = 3.5,aes(x= Network_property,Varused,fill=Network_prop))+
  geom_bar(stat ="identity")+ theme(text = element_text(size = 16),axis.text = element_text(face="bold"),
                                    legend.key.size = unit(1.2, 'cm'), #change legend key size
                                    legend.key.height = unit(1.0, 'cm'), #change legend key height
                                    legend.key.width = unit(1.2, 'cm'), #change legend key width
                                    legend.title = element_text(size=20),legend.text = element_text(size=20))
q + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#actual plot just above this#
#we see that avgclose, avghubscore and avgeigen ecn is 0 so we are removing these network properties from the training
#we realise that these features have 0 value on test and not in train data, so we dont remove
tr = diff.bal[,-c(3,7,8)]
te = valdt2[,-c(3,7,8)]
tr = diff.bal[,-c(3,7,8,9,10,11,12,14,16,17,1)]
te = valdt2[,-c(3,7,8,9,10,11,12,14,16,17,1)]
tr = diff.bal[,c(13,18)]
te = valdt2[,c(13,18)]
mod.sub = randomForest(gi ~ . , data = tr)
mod.sub = randomForest(gi ~ . , data = tr, importance = T)
pred = predict(mod.sub, newdata = te[-15])
#pred = predict(mod.sub, newdata = te[-7])
pred_gi = predict(mod.sub, newdata = te)
te$pred_gi = pred_gi
cnfmat = table(te$gi, te$pred_gi)
cnfmat
acc = sum(diag(cnfmat)/sum(cnfmat))
#refer this video for future imp
#https://www.youtube.com/watch?v=O8Un6lQlnB8
#refer this video for coreness
#https://www.youtube.com/watch?v=8sNZ5d8eNC8
varImpPlot(mod.sub, sort = T, n.var = 17, main = "Discriminatory variables")
importance(mod.sub)
varUsed(mod.sub)
importanceOrder=order(-model1$importance)
write.csv(te, "te.avgclose.cons.n4.ecc.csv", row.names = T, quote = F)
a3 = predict(mod.sub, te, type = "vote", norm.vote = T)
write.csv(a3, "pred.avgclose.cons.n4.ecc.csv", row.names = T, quote = F)
write.csv(valdt, "valdt.csv", row.names = T, quote = F)
############################################################################
#plot to show the how often the variables for mod.sub have been used
Varused <- c(46319, 70327, 38331, 40780, 42088, 31611, 37233, 37225, 46266,
70803, 69257, 68894, 70520, 70906, 34988, 32772, 34301)
Network_prop <- c(rep("avgdeg" , 1) ,rep("avgbet" , 1) , rep("avgclos" , 1) , rep("avgcor" , 1), rep("avgcon", 1),
                  rep("avgecc" , 1), rep("avgeigcen" , 1), rep("avghub" , 1), rep("avgn1" , 1), rep("avgn2" , 1),
                  rep("avgn3" , 1), rep("avgn4" , 1), rep("avgn5" , 1), rep("avgn6" , 1), rep("shortpath" , 1),
                  rep("cohesn" , 1), rep("adhesn" , 1))
data <- data.frame(Varused, Network_prop)    
Network_property = reorder(Network_prop,-Varused)
ggplot(data,add = T, las = 1,cex.axis = 1.8,aes(x= Network_property,Varused,fill=Network_prop))+
  geom_bar(stat ="identity")+ theme(legend.key.size = unit(0.8, 'cm'), #change legend key size
                                    legend.key.height = unit(0.6, 'cm'), #change legend key height
                                    legend.key.width = unit(0.6, 'cm'), #change legend key width
                                    legend.title = element_text(size=20),legend.text = element_text(size=15))
#############################################################################
#we plot the box plot of each feature to better understand the features
#avgdeg
train_nb = diff.bal
GI = reorder(train_nb$gi,-train_nb$avgdeg,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgdeg, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgdeg, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgbet
GI = reorder(train_nb$gi,-train_nb$avgbet,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgbet, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgbet, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgneighborhood5
GI = reorder(train_nb$gi,-train_nb$avgneighbor5,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor5, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor5, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgneighborhood4
GI = reorder(train_nb$gi,-train_nb$avgneighbor4,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor4, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor4, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgneighborhood1
GI = reorder(train_nb$gi,-train_nb$avgneighbor1,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor1, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor1, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgneighborhood2
GI = reorder(train_nb$gi,-train_nb$avgneighbor2,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor2, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor2, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgneighborhood3
GI = reorder(train_nb$gi,-train_nb$avgneighbor3,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor3, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor3, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgclose
GI = reorder(train_nb$gi,-train_nb$avgclose,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgclose, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgclose, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgcons
GI = reorder(train_nb$gi,-train_nb$avgcons,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgcons, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgcons, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgcore
GI = reorder(train_nb$gi,-train_nb$avgcore,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgcore, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgcore, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgeigen_cen
GI = reorder(train_nb$gi,-train_nb$avgeigen_cen,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgeigen_cen, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgeigen_cen, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avghubscore
GI = reorder(train_nb$gi,-train_nb$avghub_score,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avghub_score, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avghub_score, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
GI = reorder(train_nb$gi,-train_nb$avgneighbor6,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor6, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor6, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#avgecc
GI = reorder(train_nb$gi,-train_nb$avgecc,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgecc, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$avgecc, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#shortest_path
GI = reorder(train_nb$gi,-train_nb$shortest_path,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=shortest_path, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$shortest_path, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#adhesion
GI = reorder(train_nb$gi,-train_nb$adhesion,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=adhesion, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$adhesion, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#cohesion
GI = reorder(train_nb$gi,-train_nb$cohesion,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=cohesion, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  scale_y_continuous(limits = quantile(train_nb$cohesion, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))
#tryin calculation of network properties with a dummy network
net = read.csv("/home/nikola/brain-ppi-revise/net.csv")
net = read.csv("/home/nikola/brain-ppi-revise/net2.csv")
library(igraph)
Net <- graph_from_data_frame(net, directed=FALSE)
net = as.data.frame(degree(Net))
net$names <- rownames(net)
net$degree = degree(Net)
net = net[, -1]
row.names(net) = NULL
net$betweenness = betweenness(Net)
net$coreness = coreness(Net)
net$neighbor1 <- neighborhood.size(Net, order = 1)
net$neighbor2 <- neighborhood.size(Net, order = 2)
#matching the row name sto get the gene ids
a1 = read.csv("valdt.csv")
a2 = read.csv("te.avgclose.cons.n4.ecc.csv")
a3 = read.csv("pred.avgclose.cons.n4.ecc.csv")
b1 = merge(a1, a2, by='X')
b2 = merge(b1, a3, by = 'X')
b3 = b2[, c(2,3,6,8,9,15,26,27,28,29,30)]
#donut plot to represent the number of SL, SV, NOT
#https://r-graph-gallery.com/128-ring-or-donut-plot.html
data <- data.frame(
  category=c("SL","NOT","SV"),
  count=c(4256,6852, 8331)
)
data <- data.frame(
  category=c("SL","NOT","SV"),
  count=c(15451,24657, 2777)
)
# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  scale_fill_manual(values=c("#619CFF", "#F8766D", "#00BA38")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
#using wilcoxon test to select the features out of avgbet, deg, core, avgn1,n5
#to check a significant difference betwen SL and SV pairs
#wilcoxon test
#avgdeg
xyz3 = diff.bal
sl = xyz3[which(xyz3$gi == "SL"),]
sv = xyz3[which(xyz3$gi == "SV"),]
sl[4257:8331,] = ""
sl2 = as.numeric(sl[,1])
sv2 = as.numeric(sv[,1])
wilcox_test = wilcox.test(sl2, sv2, paired = TRUE)
#V = 4682986, p-value = 0.001095
#avgbet
sl = xyz3[which(xyz3$gi == "SL"),]
sv = xyz3[which(xyz3$gi == "SV"),]
sl[4257:8331,] = ""
sl3 = as.numeric(sl[,2])
sv3 = as.numeric(sv[,2])
wilcox_test2 = wilcox.test(sl3, sv3, paired = TRUE)
#V = 4251800, p-value = 0.0007056
#avgcore
sl = xyz3[which(xyz3$gi == "SL"),]
sv = xyz3[which(xyz3$gi == "SV"),]
sl[4257:8331,] = ""
sl4 = as.numeric(sl[,4])
sv4 = as.numeric(sv[,4])
wilcox_test3 = wilcox.test(sl4, sv4, paired = TRUE)
#V = 4441614, p-value = 0.009872
#avgneighbor1
sl = xyz3[which(xyz3$gi == "SL"),]
sv = xyz3[which(xyz3$gi == "SV"),]
sl[4257:8331,] = ""
sl5 = as.numeric(sl[,9])
sv5 = as.numeric(sv[,9])
wilcox_test4 = wilcox.test(sl4, sv4, paired = TRUE)
#V = 4441614, p-value = 0.009872
#avgneighbor5
sl = xyz3[which(xyz3$gi == "SL"),]
sv = xyz3[which(xyz3$gi == "SV"),]
sl[4257:8331,] = ""
sl6 = as.numeric(sl[,13])
sv6 = as.numeric(sv[,13])
wilcox_test5 = wilcox.test(sl5, sv5, paired = TRUE)
#V = 4677943, p-value = 0.00166
#trying to plot roc curve
library(multiROC)
set.seed(123456)
main_fltr2 = diff.bal
main_fltr2[["gi"]] = factor(main_fltr2[["gi"]])
total_number <- nrow(main_fltr2)
train_idx <- sample(total_number, round(total_number*0.6))
train_df <- main_fltr2[train_idx, ]
test_df <- main_fltr2[-train_idx, ]

rf_res <- randomForest::randomForest(gi~., data = train_df, ntree = 100)
rf_pred <- predict(rf_res, test_df, type = 'prob') 
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")

mn_res <- nnet::multinom(gi ~., data = train_df)
mn_pred <- predict(mn_res, test_df, type = 'prob')
mn_pred <- data.frame(mn_pred)
colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")

library(dummies)
true_label <- dummies::dummy(test_df$gi, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)
#final_df2 = final_df[, c(1,4,7,2,5,8,3,6,9)]
roc_res <- multi_roc(final_df, force_diag=F)
pr_res <- multi_pr(final_df, force_diag=F)

plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)

require(ggplot2)
ggplot(plot_roc_df, print.auc = T, print.auc.x = 45, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        #legend.justification=c(1, 0), legend.position=c(.99, .08),
        legend.title=element_blank(),
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
#avgneighborhood5
GI = reorder(train_nb$gi,train_nb$avgneighbor5,na.rm = TRUE)
train_nb %>%
  ggplot(aes(x=GI, y=avgneighbor5, fill=gi)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF")) +
  scale_y_continuous(limits = quantile(train_nb$avgneighbor5, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))

ggplot(train_nb, aes(x = gi, y = avgneighbor5, fill = GI)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                    labels = c("SL", "SV", "NOT")) +
  xlab("") +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(fill=guide_legend(reverse=F)) +
  scale_x_discrete(limits=c("SL", "SV", "NOT"))

library(dplyr)
library(forcats)
train_nb = file
GI = reorder(train_nb$gi,-train_nb$avgneighbor5,na.rm = TRUE)
train_nb %>% arrange(avgneighbor5) %>%
  mutate(name = factor(GI, levels=c("SL", "SV", "NOT"))) %>%
  ggplot(aes(x=gi, y=avgneighbor5, fill=gi)) + 
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")+
  scale_y_continuous(limits = quantile(train_nb$avgneighbor5, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))+guides(fill=guide_legend(reverse=F)) 
  #scale_x_discrete(limits=c("SL", "SV", "NOT")) 

GI = reorder(train_nb$gi,-train_nb$avgneighbor5,na.rm = TRUE)
train_nb %>%
  class <- factor(c("b", "b", "a", "c", "c", "c")) %>%
  ggplot( aes(x=class, y=avgneighbor5, fill=GI)) + 
  geom_boxplot() +
  xlab("class") +
  theme(legend.position="none") +
  xlab("") +
  xlab("")

file1 = train_nb[which(train_nb$gi == "SL"),]
file2 = train_nb[which(train_nb$gi == "SV"),]
file3 = train_nb[which(train_nb$gi == "NOT"),]
file = rbind(file1,file2,file3)
file %>% 
  ggplot(aes(x= fct_inorder(gi,avgneighbor5), y=avgneighbor5, fill=gi)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")
  #geom_jitter(width=0.1,alpha=0.2)  +
 

#GI = reorder(train_nb$gi,train_nb$avgneighbor5,na.rm = TRUE)
#final plot
#avgneighbor5
train_nb$gi = factor(train_nb$gi, levels = c("SL", "SV", "NOT"))
f=train_nb$gi
train_nb %>% 
     ggplot(aes(x= fct_infreq(f), y=avgneighbor5, fill=gi)) +
     geom_boxplot(outlier.shape = NA) +
     scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")+
  scale_y_continuous(limits = quantile(train_nb$avgneighbor5, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))+
  scale_x_discrete(limits=c("SL", "SV", "NOT")) 
#avgcore
train_nb %>% 
  ggplot(aes(x= fct_infreq(f), y=avgcore, fill=gi)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")+
  scale_y_continuous(limits = quantile(train_nb$avgcore, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))+
  scale_x_discrete(limits=c("SL", "SV", "NOT")) 
#avgdeg
train_nb %>% 
  ggplot(aes(x= fct_infreq(f), y=avgdeg, fill=gi)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")+
  scale_y_continuous(limits = quantile(train_nb$avgdeg, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))+
  scale_x_discrete(limits=c("SL", "SV", "NOT")) 
#avgbet
train_nb %>% 
  ggplot(aes(x= fct_inorder(f), y=avgbet, fill=gi)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF"), labels = c("SL", "SV", "NOT")) + xlab("gi")+
  scale_y_continuous(limits = quantile(train_nb$avgbet, c(0.1, 0.9))) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=20),legend.text = element_text(size=15))+
  scale_x_discrete(limits=c("SL", "SV", "NOT")) 
