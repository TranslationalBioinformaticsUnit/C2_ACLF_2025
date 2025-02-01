# Load libraries

# A) Cohort study description: Created using Adobe Illustrator.



# B)
 a <- ggplot(aclara_signature, aes(x=CLIFAD, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "CLIF-C AD", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="none")

b <- ggplot(predict_signature, aes(x=MELD, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c("gold1", "darkorange1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "MELD", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="none")

c <- ggplot(predict_signature, aes(x=CLIFAD, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c("gold1", "darkorange1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "CLIF-C AD", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="none")

d <- ggplot(aclara_signature, aes(x=MELD, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "MELD", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="none")

# C)

a <- ggplot(data = predict_signature, aes(x=gene, y=values)) + geom_boxplot(aes(fill=groups))+
  scale_fill_manual(values=c( "gold1", "darkorange1", "firebrick4")) + labs(x = "", y = "C2 Signature", title = "PREDICT")+
  facet_wrap( ~ gene, scales="free") + theme(legend.position="none")

b <- ggplot(data = aclara_signature, aes(x=gene, y=values)) + geom_boxplot(aes(fill=groups))+
  scale_fill_manual(values=c( "gold1", "firebrick4")) + labs(x = "", y = "C2 Signature", title = "ACLARA")+
  facet_wrap( ~ gene, scales="free") + theme(legend.position="none")


# D)

a <-  ggplot(data = predict_signature, aes(x=Mortality, y=GeneSet)) + geom_boxplot(aes(fill=Mortality))+
  labs(x = "", y = "C2 signature", title = "", fill="") +
  theme(axis.text=element_text(size=10)))
  )

c <-  ggplot(data = aclara_signature, aes(x=Mortality, y=GeneSet)) + geom_boxplot(aes(fill=Mortality))+
  labs(x = "", y = "C2 signature", title = "", fill="") +
  theme(axis.text=element_text(size=10)))
  )

# E)

plot.roc(roc_obj_predict_A,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC ACLF (28d)" )
lines(roc_obj_aclara_A, col = "#0072B2")
legend("bottomright", legend = c("PREDICT (AUC=0.657)", "ACLARA (AUC=0.637)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)


plot.roc(roc_obj_predict_M,print.auc=F,
         col="#D55E00",xlab="1-Specificity",ylab="Sensitivity", main ="ROC Mortality (90d)" )
lines(roc_obj_aclara_M, col = "#0072B2")
legend("bottomright", legend = c("PREDICT (AUC=0.659)", "ACLARA (AUC=0.625)"), col = c("#D55E00", "#0072B2"),  lwd = 2, box.col  = "black", xjust = 1, yjust = 2,cex=0.8)
