# load libraries
library(pROC)


plot.roc(roc_obj_predict_A,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC ACLF (28d)" )
lines(roc_obj_predict_B, col = "#0072B2")
legend("bottomright", legend = c("CLIF-C AD (AUC=0.657)", "CLIF-C AD + C2 (AUC=0.637)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)


plot.roc(roc_obj_predict_C,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC Mortality (90d)" )
lines(roc_obj_predict_D, col = "#0072B2")
legend("bottomright", legend = c("CLIF-C AD (AUC=0.659)", "CLIF-C AD + C2 (AUC=0.625)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)




plot.roc(roc_obj_aclara_A,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC ACLF (28d)" )
lines(roc_obj_aclara_B, col = "#0072B2")
legend("bottomright", legend = c("CLIF-C AD (AUC=0.657)", "CLIF-C AD + C2 (AUC=0.637)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)


plot.roc(roc_obj_aclara_C,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC Mortality (90d)" )
lines(roc_obj_aclara_D, col = "#0072B2")
legend("bottomright", legend = c("CLIF-C AD (AUC=0.659)", "CLIF-C AD + C2 (AUC=0.625)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)

