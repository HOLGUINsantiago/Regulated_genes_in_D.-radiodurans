
library(dplyr)
library(FactoMineR)
library(factoextra)


df <- read.csv("D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/Compte rendu/pv_par_gene.csv.csv")
associated_genes452 <- read.delim("D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/Compte rendu/ChipSeq_motifs/associated_genes452.tsv")
associated_genes136 <- read.delim("D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/Compte rendu/ChipSeq_motifs/associated_genes136.tsv")

genes_buscar <- c("DRO_A0153", "DRO_A0342", "DRO_0003", "DRO_0070", "DRO_0099",
                  "DRO_0171", "DRO_0219", "DRO_0323", "DRO_0421", "DRO_0596",
                  "DRO_0657", "DRO_0899", "DRO_1033", "DRO_1140", "DRO_1255",
                  "DRO_1280", "DRO_1673", "DRO_1751", "DRO_1755", "DRO_1891",
                  "DRO_1899", "DRO_2230", "DRO_2249", "DRO_2308")


df_associated <- associated_genes452 %>%
  select(Gene_Name, Region) %>%
  rename(Gene = Gene_Name)


df_summary <- df %>%
  group_by(Gene) %>%
  summarise(
    FC_sup_2.5_D37 = sum(log2FC > 2.5 & Echantillon == "D37"),
    FC_sup_1.5_D37 = sum(log2FC > 1.5 & Echantillon == "D37"),
    FC_sup_2_4fois_D37 = sum(log2FC > 2 & Echantillon == "D37") >= 4,
    
    FC_inf_1.5_W37 = sum(log2FC < 1.5 & Echantillon == "W37"),
    FC_inf_1_4fois_W37 = sum(log2FC < 1 & Echantillon == "W37") >= 4,
    Significatif_diff = sum(P_adj < 0.01 & Echantillon == "D37" & log2FC > 1) - sum(P_adj < 0.01 & Echantillon == "W37" & log2FC > 1),
    Associated = Gene %in% associated_genes452$Gene_Name,
    Reference = Gene %in% genes_buscar
  ) %>%
  left_join(df_associated, by = "Gene") %>%  
  mutate(Associated_region = ifelse(is.na(Region), "NOT", Region)) %>%  
  select(-Region)%>%  
  distinct(Gene, .keep_all = TRUE)


df_filtered <- df_summary %>%
  mutate(
    Condition = case_when(
      FC_sup_2.5_D37 >= 1 & FC_inf_1_4fois_W37 & FC_inf_1.5_W37 ==5 & Significatif_diff >= 3  ~ 1,  
      FC_sup_1.5_D37 >= 3 & FC_sup_2.5_D37 >= 2 & FC_inf_1.5_W37 >=3 & Significatif_diff >= 3 ~ 2,
      FC_sup_1.5_D37 >= 1 & FC_inf_1_4fois_W37  & Associated & Significatif_diff >= 1 ~ 3,  
      FC_sup_2_4fois_D37 & FC_inf_1_4fois_W37  & Significatif_diff >= 3 ~ 4,  
      TRUE ~ NA_real_  
    )
  ) %>%
  filter(!is.na(Condition)) %>%  
  select(Gene, Condition, Reference) %>%
  distinct(Gene, .keep_all = TRUE)


genes_L2 <- c("DRO_0595", "DRO_0614", "DRO_1666", "DRO_1769", "DRO_2182", "DRO_2578")
genes_L3 <- c("DRO_A0343", "DRO_C0016", "DRO_0004", "DRO_0216", "DRO_0422", "DRO_0558",
              "DRO_0898", "DRO_1032", "DRO_1881", "DRO_1890", "DRO_2043", "DRO_2307", "DRO_2414")
genes_L4 <- c("DRO_A0058", "DRO_A0214", "DRO_A0326", "DRO_B0035", "DRO_C0011", "DRO_C0033",
              "DRO_0199", "DRO_1280", "DRO_1459", "DRO_1547", "DRO_1673", "DRO_1685",
              "DRO_1829", "DRO_1899", "DRO_2399", "DRO_2429")


df_summary <- df_summary %>%
  left_join(df_filtered %>% select(Gene, Condition), by = "Gene") %>%
  mutate(
    Passe_Filtre = !is.na(Condition),
    Condition = ifelse(is.na(Condition), -2, Condition),
    L2 = as.numeric(Gene %in% genes_L2),
    L3 = as.numeric(Gene %in% genes_L3),
    L4 = as.numeric(Gene %in% genes_L4),
    article=as.numeric(L2 ==1 | L3 ==1 | L4 ==1)
  )


genes_non_filtres <- df_summary %>%
  filter((L2 ==1 | L3 ==1 | L4 ==1) & Passe_Filtre == FALSE)


write.csv(df_summary, "D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/Compte rendu/resultats_filtres.csv", row.names = FALSE)
write.csv(df_filtered, "D:/Recherche/Master/M1 - Neurosciences - Physiologie Paris-Saclay/Bloc mineur/ADMO/Compte rendu/genes_filtres.csv", row.names = FALSE)



df_numeric <- df_summary %>%
  ungroup() %>%  
  mutate(across(where(is.logical), as.numeric)) %>%
  select(-Gene) %>%  
  select(where(is.numeric))  



pca_res <- PCA(df_numeric, scale.unit = TRUE, graph = FALSE)


fviz_pca_biplot(pca_res, label = "var", habillage = df_summary$Passe_Filtre,
                addEllipses = TRUE, ellipse.level = 0.95, repel = TRUE)


fviz_pca_biplot(pca_res, label = "var", habillage = df_summary$article,
                addEllipses = TRUE, ellipse.level = 0.95, repel = TRUE)

df_numeric_clean <- df_summary %>%
  ungroup() %>% 
  select(-c(Gene, Condition, L2, L3, L4, article, Passe_Filtre)) %>%
  mutate(across(where(is.logical), as.numeric)) %>%
  select(where(is.numeric))  

pca_res_clean <- PCA(df_numeric_clean, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(pca_res_clean, label = "var", 
                addEllipses = TRUE, ellipse.level = 0.95, repel = TRUE)

coord_ind_clean <- as.data.frame(pca_res_clean$ind$coord)
rownames(coord_ind_clean) <- df_summary$Gene  


ref_coord_clean <- coord_ind_clean[df_summary$Reference == TRUE, ]
filtre_coord_clean <- coord_ind_clean[df_summary$Passe_Filtre == TRUE, ]
l234_coord_clean <- coord_ind_clean[df_summary$L2 == 1 | df_summary$L3 == 1 | df_summary$L4 == 1, ]

distance_ref_filtre_clean <- mean_distance(ref_coord_clean, filtre_coord_clean)
distance_ref_l234_clean <- mean_distance(ref_coord_clean, l234_coord_clean)

print(paste("Distance Reference - Passe_Filtre (sans biais) :", distance_ref_filtre_clean))
print(paste("Distance Reference - L2-4 (sans biais) :", distance_ref_l234_clean))

