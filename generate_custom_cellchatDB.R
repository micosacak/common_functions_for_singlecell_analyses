convert2csv = function(the_database){
  #the_database <- the_database.zebrafish # set the_database <- the_database.human if working on the human dataset
  interaction_input <- the_database$interaction
  complex_input <- the_database$complex
  cofactor_input <- the_database$cofactor
  geneInfo <- the_database$geneInfo
  write.csv(interaction_input, file = "interaction_input_the_database.csv")
  write.csv(complex_input, file = "complex_input_the_database.csv")
  write.csv(cofactor_input, file = "cofactor_input_the_database.csv")
  write.csv(geneInfo, file = "geneInfo_input_the_database.csv")
}
convert2db = function(){
  options(stringsAsFactors = FALSE)
  interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
  complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
  cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
  geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
  new_database <- list()
  new_database$interaction <- interaction_input
  new_database$complex <- complex_input
  new_database$cofactor <- cofactor_input
  new_database$geneInfo <- geneInfo
  return(new_database)
}
