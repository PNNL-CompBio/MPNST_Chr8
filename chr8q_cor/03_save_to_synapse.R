##Save all analysis and data results to synapse


new_syn_upload <-function(file_list, synId){
  for(i in file_list)
    if(isDirectory(i)){
      res <- synapser::synStore(synapser::Folder(i, parentId=synId))
      print(paste('storing directory',i))
      setwd(i)
      res <- new_syn_upload(list.files('.'), res$id)
      setwd('..')
    }else{
      res <- synapser::synStore(synapser::File(i,parentId=synId))
      print(paste('storing file',i))
    }
  return(res)
}



setwd(base.path)
setwd('..')

top.syn <- "syn60219614"

#SG replaced uplod with simpler recursion
new_syn_upload(file_list=basename(base.path),synId=top.syn)