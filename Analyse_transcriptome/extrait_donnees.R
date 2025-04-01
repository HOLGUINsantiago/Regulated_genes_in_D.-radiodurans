folders <- c("D37", "w37")

base_path <- "D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/multi res/"

base_path <- getwd()
base_path <- paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")

process_files <- function(folder) {
  files <- list.files(paste0(base_path, folder, "/"), pattern="*.tabular", full.names=TRUE)
  
  df_list <- lapply(files, function(f) {
    df <- read.delim(f, header=FALSE)
    colnames(df) <- c("Gene", "Base mean", "log2FC", "StdErr", "Wald-Stats", "P-value", "P_adj")
    
    temps <- tools::file_path_sans_ext(basename(f))
    
    df$Temps <- temps
    df$Echantillon <- folder
    
    return(df)
  })
  
  return(do.call(rbind, df_list))
}

df_combine <- do.call(rbind, lapply(folders, process_files))

df_def <- df_combine[,c( "Gene","log2FC", "P_adj", "Temps", "Echantillon")]

write.csv(df_def, paste0(base_path, "pv_par_gene.csv"), row.names=FALSE)


