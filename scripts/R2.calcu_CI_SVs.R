#!R
## Rscriotp to plot bp precision distribution and calculate CI of callers
args = commandArgs(trailingOnly=TRUE)

ppre=paste(args[1],'/',sep='')
data_file=paste(ppre,'STEP1_bp_off_ILL_vs_PB_both_end_separate_RC50.txt',sep='')

data=read.table(data_file)
#par(mfrow=c(3,5))
out_table=data.frame('caller'=0,'left'=0,'right'=0)
rec=0
for (k1 in unique(data[,1])){
	png(paste(paste(ppre,k1,sep=''),'consensus_bp_off_illumina_vs_pacbio_both_end_separate_based_on_RC50','png',sep='.'))
	temp=data[data[,1]==k1,]
	left_quantile=quantile(temp[,3], c(.05,.95)) 
	right_quantile=quantile(temp[,4], c(.05,.95)) 
	left_keep_quantile=quantile(temp[,3], c(.1,.9)) 
	right_keep_quantile=quantile(temp[,4], c(.1,.9)) 
	rec=rec+1
	out_table[rec,]=c(k1,range(c(left_keep_quantile,right_keep_quantile)))
	xrange=range(c(left_quantile,right_quantile))
	hist(temp[temp[,3]>left_quantile[1]&temp[,3]<left_quantile[2],][,3],cex.axis=1.5,xlab='',ylab='',cex.axis=2,xlim=xrange,col='blue',breaks=100,main=k1,cex.main=2)
	hist(temp[temp[,4]>right_quantile[1]&temp[,4]<right_quantile[2],][,4],xlab='',ylab='',col='red',cex.axis=2,breaks=100,add=T)
	dev.off()
	}
	png(paste('legend','consensus_bp_off_illumina_vs_pacbio_both_end_separate_based_on_RC50','png',sep='.'))
	plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n')
	legend("topleft",cex=2, c("left breakpoints", "right breakpoints"), fill=c("blue", "red"),bty='n')	
	dev.off()

write.table(out_table, paste(ppre,'consensus_bp_CI_90Quantile.txt',sep=''), quote=F, col.names=F, row.names=F, sep='\t')

