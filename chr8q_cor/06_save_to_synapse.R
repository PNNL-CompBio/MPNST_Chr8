##Save all analysis and data results to synapse

library(synapser)
synLogin()
new_syn_upload <-function(file_list, synId){
  for(i in file_list)
    if(file_test('-d',i)){
      if(length(grep('\\.',i)!=0)){ ##these are incorrect
        print(paste('bad directory',i))
        file.remove(i)
        next
      }
      res <- synapser::synStore(synapser::Folder(i, parentId=synId))
      print(paste('storing directory',i))
      setwd(i)
      if(length(list.files('.'))>0)
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