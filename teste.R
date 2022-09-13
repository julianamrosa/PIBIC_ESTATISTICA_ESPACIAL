ggplot(out,aes(x=as.numeric(out$h1),y=as.numeric(out$res1))) + geom_point(colour="red",alpha=0.5,size=3) +
  scale_y_continuous() + 
  scale_x_continuous() + 
  geom_point(aes(x=as.numeric(out$h2),y=as.numeric(out$res2)),alpha=.5,size=3, colour="darkblue") + 
  labs(x= "Bandwidth", y="aic") +
  theme_test() +
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

q1<- quantile(c(as.numeric(out$res1), as.numeric(out$res2)), 0.25, na.rm=T)
q3 <- quantile(c(as.numeric(out$res1), as.numeric(out$res2)), 0.75, na.rm=T)
iqr <- q3-q1

par(mgp = c(1.5, 0.4, 0), tcl = -0.25, mar=c(3, 2.4, 1.3, 0.6))
plot(x=as.numeric(out$h1), y=if(sum(c(as.numeric(out$res1[!is.na(out$res1)]),
                                      as.numeric(out$res2[!is.na(out$res2)]))<=q1-1.5*iqr
                                    | c(as.numeric(out$res1[!is.na(out$res1)]),
                                        as.numeric(out$res2[!is.na(out$res2)]))>=q3+1.5*iqr)>0){log10(as.numeric(out$res1)); logaritmo <- TRUE}else{as.numeric(out$res1) logaritmo <- FALSE},
type="p", pch=ifelse(out$res1[!is.na(out$res1)]==min(c(out$res1[!is.na(out$res1)],
                                                       out$res2[!is.na(out$res2)])), 17, 19),
yaxt="n", xaxt = "n", xlab="Bandwidth", ylab="aic",
     cex.main=1, cex.lab=0.9, cex.axis=0.7,
     col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5),
     cex=1.3, alpha=0.5)
points(x=as.numeric(out$h2), y=if(sum(c(as.numeric(out$res1[!is.na(out$res1)]),
                                        as.numeric(out$res2[!is.na(out$res2)]))<=q1-1.5*iqr
                                      | c(as.numeric(out$res1[!is.na(out$res1)]),
                                          as.numeric(out$res2[!is.na(out$res2)]))>=q3+1.5*iqr)>0){log10(as.numeric(out$res2))}else{as.numeric(out$res2)}, type="p", pch=ifelse(out$res2[!is.na(out$res2)]==min(c(out$res1[!is.na(out$res1)], out$res1[!is.na(out$res2)])), 17, 19),
       xaxt = "n", col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5), cex=1.3, alpha=0.5)
axis(1, cex.axis=0.8)
axis(side = 2, lwd = 0, lwd.ticks = 2, las = 2, cex.axis=0.8)

