library(ggplot2)

#input_data=input_data[1:16,]
for (group in 1:40){
  input_file_name=paste("C:/Research_data/mtDNA_identification/exp1_group_result/all_sample_mtDNA_summary/group",group,"all_sample_summary.txt",sep="")
  input_data=read.table(input_file_name,header=TRUE)
  output_png_name=paste("C:/Research_data/mtDNA_identification/exp1_group_result/all_sample_mtDNA_png/group",group,"all_sample_summary.png",sep="")
  png(output_png_name,height = 632,width = 1231)
  p=ggplot(data = input_data,mapping=aes(x = SampleName,y= Count,fill=Category))
  p=p+geom_bar(stat="identity")
  p=p+theme_light()+theme(axis.text.x = element_text(angle = 90,vjust=0.4))+ggtitle(paste("Group",group,"Reads Count Per Sample",sep=" "))+ylab("Read Counts")
  p=p+ylim(0,81000)
  print(p)
  dev.off()
}
