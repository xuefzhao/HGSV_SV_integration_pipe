#!R
##Rscript to calcuate confident intervals of insertion points

plot_and_quantile<-function(k1,ppre,fileout){
	png(paste(ppre, paste('Ins_dis',k1,'png',sep='.'),sep=''))
	data=read.table(paste(ppre, paste('Ins_dis',k1,'rec',sep='.'),sep=''))
	data=data[abs(data[,1])<1000,]
	quantile=quantile(data[,1], c(.1,.9)) 
	par(mfrow=c(2,1))
	x=data[,1]
	hist.data = hist(x, plot=F,breaks=100)
	hist.data$counts = log10(hist.data$counts)
	plot(hist.data,ylim=c(0,max(hist.data$counts )),main=k1,xlab='distance of insert points')
	x=data[,2]
	hist.data = hist(x, plot=F,breaks=100)
	hist.data$counts = log10(hist.data$counts)
	plot(hist.data,ylim=c(0,max(hist.data$counts )),main='',xlab='difference between insert length')
	dev.off()	
	return(c(k1,quantile))
	}

caller_name_readin<-function(ppre){
	caller_name=c()
	for(i in list.files(ppre)){
		if(strsplit(i,'[.]')[[1]][1]=='Ins_dis' & strsplit(i,'[.]')[[1]][3]=='rec'){
			caller_name=c(caller_name, strsplit(i,'[.]')[[1]][2])
		}
	}
	return(caller_name)
	}

args = commandArgs(trailingOnly=TRUE)

ppre=paste(args[1],'/',sep='')
caller_name=caller_name_readin(ppre)
out_frame=data.frame('algorithm'=0,'quantile_10'=0,'quantile_90'=0)
rec=0

for(k1 in caller_name){	
	rec=rec+1
	test=plot_and_quantile(k1,ppre,fileout)	
	out_frame[rec,]=test
	}	

write.table(out_frame, paste(ppre,'consensus_INS_CI_90Quantile.txt',sep=''), quote=F, sep='\t', col.names=F, row.names=F)


