#Test gprofiler2_analysis function with merge_result and clusterProfiler::dotplot visualisation function

#Creating the list of named date.frames with your data
des_res <- list(
   df1 = read_csv(system.file("extdata/res_5ari_signif_down", package = "rnaseq")),
   df2 = read_csv(system.file("extdata/res_diet_signif_down", package = "rnaseq")),
   df3 = read_csv(system.file("extdata/res_5ari_signif_up", package = "rnaseq")),
   df4 = read_csv(system.file("extdata/res_diet_signif_up", package = "rnaseq"))
   )

#Generatin the enrichment result objects that will be contained in enr$gp and enr$cp
enr <- gprofiler2_analysis(de_res)

#Creating a list of the enrichResult objects we want to visualise
enr_list <- list(enr$cp$df1$enrichResult,
                 enr$cp$df2$enrichResult,
                 enr$cp$df3$enrichResult,
                 enr$cp$df4$enrichResult)

#Name each object for a corresponding name in the x axis of the dotplot
names(enr_list) <- c("enrari_down",
                     "enrdiet_down",
                     "enrari_up",
                     "enrdiet_up")

#Merge the enrichResult object for visualisation with merge_result
merged_result <- clusterProfiler::merge_result(enr_list)

#Generate a dotplot from the merge_result object
clusterProfiler::dotplot(merged_result)
