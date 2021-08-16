suppressMessages(library(stringr))
suppressMessages(library(Biostrings))


get_closest_motif<-function(x){
	transid=bed[x,1]
	pos=as.integer(bed[x,3])
	seq=trans_seq[[transid,exact = FALSE]]
	if(length(seq)==0){
		return(c(NA,NA,NA))
	}	
	base=as.character(as.data.frame(subseq(seq, start=pos,end=pos))[1,1])
	tmp=matchPattern("RRACH", seq, fixed = FALSE)# mapping the motif to seuqence of transcript
	tmp2=as.vector(tmp@ranges@start+2) #obtain position of center A
	if(length(tmp2)==0){
		return(c(NA,NA,NA))
	}
	index=which.min(abs(as.numeric(pos)-tmp2)) #get index of closest position
	return(c(tmp2[index],as.data.frame(tmp)[index,],base)) #return the closest position and instantiated motif (like 'AAACA')
}

args=commandArgs(T)
file_path=args[1]
print(file_path)
bed=read.delim(file_path,header=F)
motif=args[2]
fasta=args[3]
out_preifx=args[4]
out_path=args[5]

trans_seq <- readDNAStringSet(fasta, "fasta")

parallel_process<-function() 
{
library(parallel)
clnum<-detectCores() 
mc <- getOption("mc.cores", clnum%/%2)

res <- mclapply(1:nrow(bed), get_closest_motif, mc.cores = mc)
result2=data.frame(matrix(unlist(res), ncol = 3,byrow = T),stringsAsFactors = FALSE)
names(result2)=c('m6a_pos','motif','error_base')
final_res=cbind(bed,result2)
final_res$dist2motif=final_res$V3-as.numeric(final_res$m6a_pos)
setwd(out_path)
write.table(final_res,file=paste(out_preifx,'closet_motif.txt',sep="_"),sep='\t',quote=F,col.names=F,row.names=F)
}
parallel_process()

