liste_output_complet <- list.files(path = path_output, pattern = (".tmp"))# on va traiter l'ensemble des fichiers presents dans le dossier 1.extraction
fichier_sortie <- paste(path_output,'sortie_f.txt',sep='')
if(length(liste_output_complet)!=0){
      content <- matrix(scan(paste(path_output,liste_output_complet[1],sep=''),what=character(0),skip=0,quiet=TRUE),ncol=14,byrow=TRUE)
      if(length(liste_output_complet)>=2){
            for(num_fichier in 2:length(liste_output_complet)){
                  content <- rbind(content,matrix(scan(paste(path_output,liste_output_complet[num_fichier],sep=''),what=character(0),skip=1,quiet=TRUE),ncol=14,byrow=TRUE))
            }
      }

      write.matrix(content,file = fichier_sortie,sep='\t')
}

if(system_exploitation == "windows"){
      system('cmd /c del_output_tmp.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
}
if(system_exploitation == "linux"){
      system('rm ../2.output/*.tmp')
}