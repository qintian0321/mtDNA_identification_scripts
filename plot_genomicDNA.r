library(ggplot2)

#input_data=input_data[1:16,]
for (group in c(41,42,43,46,47)){
  input_file_name=paste("C:/Research_data/mtDNA_identification/exp2_group_result/all_sample_genomicDNA_summary/group",group,"all_sample_summary.txt",sep="")
  input_data=read.table(input_file_name,header=TRUE)
  output_png_name=paste("C:/Research_data/mtDNA_identification/exp2_group_result/all_sample_genomicDNA_png/group",group,"all_sample_summary.png",sep="")
  png(output_png_name,height = 632,width = 1231)
  if (group==41){
    group_name="ETFB"
  }else if (group==42){
    group_name="DTNA"
  }else if (group==43){
    group_name="HJ"
  }else if (group==44){
    group_name="ABCC9"
  }else if (group==45){
    group_name="L1PA3"
  }else if (group==46){
    group_name="GAPDH"
  }else if (group==47){
    group_name="CYCS"
  }
  #relevel
  input_data$Chr=factor(input_data$Chr,levels=c("chr1","chr2","chr3","chr4","chr5",
                                                "chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15",
                                                "chr16","chr17","chr18","chr19","chr20",
                                                "chr21","chr22","chrX","chrY","chrM","other"))
  p=ggplot(data = input_data,mapping=aes(x = Chr,y= Count,fill=SampleName))
  p=p+geom_bar(stat="identity",position="dodge")
  p=p+theme_light()+theme(axis.text.x = element_text(angle = 90,vjust=0.4))+ggtitle(paste(group_name,"Reads Count Per Sample",sep=" "))+ylab("Read Counts")
  if (group==46 | group == 47){
    p=p+ylim(0,2000)
  }else{
    p=p+ylim(0,60000)    
  }
  print(p)
  dev.off()
}