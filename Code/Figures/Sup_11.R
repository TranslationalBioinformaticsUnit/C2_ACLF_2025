# load libraries

# A) Scatter plots

 a <- ggplot(predict_signature, aes(x=WBC, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "darkorange1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(30, 50, 70,90), limits = c(25, 95)) + 
  labs(x = "WBC", y = "C2 signature", title = "PREDICT") +
  theme_bw( ) + theme(legend.position="none")

b <- ggplot(aclara_signature, aes(x=WBC, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "WBC", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="none")

c <- ggplot(predict_signature, aes(x=CRP, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "darkorange1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(30, 50, 70,90), limits = c(25, 95)) + 
  labs(x = "CRP", y = "C2 signature", title = "PREDICT") +
  theme_bw( ) + theme(legend.position="bottom")

b <- ggplot(aclara_signature, aes(x=CRP, y=values, color=groups)) + 
  geom_point(size=0.4) + scale_color_manual(values=c( "gold1", "firebrick4")) +
  facet_wrap( ~ gene, scales="free") + theme(legend.position="bottom")+ 
  scale_x_continuous(breaks = c(25, 50, 75), limits = c(20, 80)) + 
  labs(x = "CRP", y = "C2 signature", title = " ") +
  theme_bw( ) + theme(legend.position="bottom")



# B) Box plots

a <- ggplot(data = geneSet_predict, aes(x=BINF, y=GeneSet, fill=Group)) + geom_boxplot()+
  scale_fill_manual(values=c( "gold1", "darkorange1", "firebrick4")) + 
  labs(x = " ", y = "C2 Signature", title = "PREDICT", fill="") +
  theme(legend.position="bottom",
        legend.margin = margin(t = -15, r = 0, b = 0, l = 0),
        axis.text=element_text(size=10),
        axis.text.x=element_text(size=11, vjust = 0),
        axis.text.y=element_text(size=11),plot.title = element_text(face = "bold", hjust = 0.5, size = 10)

b <- ggplot(data = geneSet_aclara, aes(x=BINF, y=GeneSet, fill=Group)) + geom_boxplot()+
  scale_fill_manual(values=c( "gold1", "firebrick4")) + 
  labs(x = "", y = "C2 Signature", title = "ACLARA", fill="") +
  theme(legend.position="bottom",
        legend.margin = margin(t = -15, r = 0, b = 0, l = 0),
        axis.text=element_text(size=10),
        axis.text.x=element_text(size=11, vjust = 0),
        axis.text.y=element_text(size=11),plot.title = element_text(face = "bold", hjust = 0.5, size = 10)
  )
