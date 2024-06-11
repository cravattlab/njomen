## Analysis of peptide proteotypicity: Fig. 3h and Extended Data Fig. 4g in Njomen et al. "Multi-Tiered Chemical Proteomic Maps of Tryptoline Acrylamide-Protein Interactions in Cancer Cells" (2024)

# Install packages if needed ---- 

    install.packages("stringr") 
    install.packages("OrgMassSpecR")
    install.packages("BiocManager")
    BiocManager::install("UniProt.ws")
    install.packages("openxlsx")

# Set working directory ----
    
    setwd('[...]/data') # replace with working directory 
    
    
# Load libraries ---- 

    library("stringr")
    library("OrgMassSpecR")
    library("UniProt.ws")
    library("openxlsx")
    
# Run analysis of proteotypicity ----
          
      # 0. Load custom functions ----
          
          # return positions of cysteines in peptide sequence (one-letter amino-acid code)
              cys.position <- function(peptide, seq.start) {

                res <- unlist(gregexpr("C", peptide)) + seq.start - 1
                res <- paste(res, collapse = ", ")

                return(res)

              }

     
      # 1. Load input data ---- 
        
            # list of liganded proteins (protein- and cysteine-directed ABPP, Ramos and 22Rv1 cells)
            proteins <- read.table(file = "liganded-proteins.txt", header = T, sep = "\t", quote = "", check.names = F) %>% as.data.frame()
              
            # list of quantified cysteine sites in liganded proteins (cysteine-directed ABPP, Ramos and 22Rv1 cells)
            sites <- read.table(file = "liganded-proteins-qsites.txt", header = T, sep = "\t", quote = "", check.names = F) %>% as.data.frame()
                sites$Accession <- paste(sites$Uniprot, sites$Site)
              
              
      # 2a. Predict proteotypicity of cysteine-containing tryptic peptides in liganded proteins ----
              
              # retrieve protein sequences from UniProt (~1 min)
                  proteins$peptide <- NA
              
                  start <- Sys.time()
                  for (i in 1:nrow(proteins)) { proteins$peptide[i] <- queryUniProt(query = paste("accession:",proteins$Uniprot[i],sep=""), fields = "sequence")} 
                  end <- Sys.time()
                  end - start
                  

              # generate table of cysteine-containing tryptic peptides
                  errors <- 0
                  sites.all <- NULL    
                  
                  for (i in 1:nrow(proteins)) {
                  
                      temp <- Digest(sequence = as.character(proteins$peptide[i]), enzyme = "trypsin", missed = 0) %>% try()
                          if ("try-error" %in% class(temp)) {
                            errors <- c(unlist(errors), i)
                            next}
                      
                      temp <- temp[grepl(pattern = "C", x = temp$peptide),]
                          if (nrow(temp) == 0) {next}
                      temp$cys <- mapply(FUN = cys.position, peptide=temp$peptide, seq.start=temp$start)
                      temp$Protein <- proteins$Protein[i]
                      temp$Uniprot <- proteins$Uniprot[i]
                      temp$length <- temp$stop - temp$start + 1
                      
                      temp <- subset(temp, select=c("Protein","Uniprot","cys","peptide", "length"))
                      
                      sites.all <- rbind(sites.all,temp)
                      
                      }
                  
                  # write input for DeepMSPeptide as txt file
                      write.table(x = sites.all$peptide, file = "deepms-input.txt", row.names = F, col.names = F, quote = F)
                  
                  
        # 2b. Add probability of prediction (DeepMSPeptide): refer to https://github.com/vsegurar/DeepMSPeptide for installation and run instructions ---- 
                  
                  # run DeepMSpeptide via terminal (MacOS/Unix) or command prompt (Windows PC): run using Python 3.9.6 (tensorflow 2.13.1, keras 2.13.1)
                      
                  # read deepMSpeptide output - ensure it is located in the R working directory
                      deep <- read.table("deepms-input_Predictions.txt", header = T) %>% as.data.frame()
                      
                  # add to previous table    
                  sites.all$detectability <- deep$Prob[match(x = sites.all$peptide, table = deep$Peptide)]
                  
                  
        # 3. Add experimental data (quantitation in cysteine-directed ABPP and hyperreactivity) ----       
                    
              # indicate sites quantitated in cysteine-directed ABPP experiments (Ramos and 22Rv1 cells)
                  sites.all$cys.abpp.detected <- F
                  
                  for (i in 1:nrow(sites.all)) { 

                    acc <- paste(sites.all$Uniprot[i], unlist(strsplit(sites.all$cys[i], split=", ")))

                    k <- which(sites$Accession %in% acc)
                    
                    if (length(k) == 0) {next} else { sites.all$cys.abpp.detected[i] <- T }
                    
                    }
                  
                  
                # add hyperreactivity data
                  
                  # read data
                  hr.file <- "20230217_CSOL_Ramos_22Rv1_8M-urea.tsv"
                      hr <- read.table(file = hr.file, header = T, sep = "\t", quote = "", check.names = F) %>% as.data.frame()

                  # for each protein entry, determine which cell line leads to higher abundance and calculate average of corresponding values under native and urea denaturation conditions
                  
                  hr$sum.ramos <- hr[,grep(pattern = "Ramos_", colnames(hr))] %>% rowSums()
                  hr$sum.22rv1 <- hr[,grep(pattern = "22Rv1_", colnames(hr))] %>% rowSums()
                  
                  hr$cell <- ""
                  hr$native <- NA
                  hr$urea <- NA
                  
                      for (i in 1:nrow(hr)) {
                        
                        if (hr$sum.ramos[i] > hr$sum.22rv1[i]) { 
                          
                          hr$cell[i] <- "Ramos"
                          hr$native[i] <- hr[i,grep(pattern = "Ramos_native", colnames(hr))] %>% unlist   %>% mean()
                          hr$urea[i] <- hr[i,grep(pattern = "Ramos_urea", colnames(hr))] %>% unlist  %>% mean()
                          
                          } else if (hr$sum.ramos[i] < hr$sum.22rv1[i]) { 
                            
                          hr$cell[i] <- "22Rv1"
                          hr$native[i] <- hr[i,grep(pattern = "22Rv1_native", colnames(hr))] %>% unlist %>% mean()
                          hr$urea[i] <- hr[i,grep(pattern = "22Rv1_urea", colnames(hr))] %>% unlist  %>% mean()  
                            
                          } else { next }
                        
                      }
                      
                  hr <- hr %>% subset(select = c(Accession, Site, Description, cell, native, urea)) 
                
                  
                  # calculate hyperreactivity
                  
                    # normalize to sum 100    
                    nf <- 100/(hr$native + hr$urea) # normalization factor
                    
                    hr$native <- hr$native*nf  
                    hr$urea <- hr$urea*nf 
                  
                    # fix smallest possible value to 5%
                    nat.adj <- lapply(X = hr$native, FUN = max, 5) %>% unlist %>% as.numeric
                    ure.adj <- lapply(X = hr$urea, FUN = max, 5) %>% unlist %>% as.numeric
                    
                    # calculate ratio and log-transform
                    hr$reac <- nat.adj/ure.adj
                    hr$reac <- -log(hr$reac, base = 2) %>% lapply(FUN = round, digits = 2) %>% unlist
                    
                  
                  # add to previous table
                  
                  sites.all$cys.denat <- ""
                  sites.all$enhancement.log2FC <- ""
                  for (i in 1:nrow(sites.all)) { 

                    acc <- paste(sites.all$Uniprot[i], unlist(strsplit(sites.all$cys[i], split=", ")))
                    k <- which(hr$Accession %in% acc)
                    
                    if (length(k) == 0) {next} else {
                      sites.all$cys.denat[i] <- hr$Accession[k] %>% lapply(FUN=word, start = 2, end = 2, sep = " ") %>% unlist %>% paste(collapse=", ")
                      sites.all$enhancement.log2FC[i] <- hr$reac[k] %>% mean } }
                  
                  
            # 4. Generate tables for manuscript figures ----
                  
                  # Figure 3h: Detection probability (DeepMSPeptide score) of tryptic peptides of stereoselectively liganded proteins (protein- or cysteine-directed ABPP)
                  
                  q <- sites.all %>% subset(cys.abpp.detected == T, select=c(Protein, Uniprot, cys, detectability))
                  nq <- sites.all %>% subset(cys.abpp.detected == F, select=c(Protein, Uniprot, cys, detectability)) 
                  
                  q.table <- matrix(data = "", ncol = 2, nrow = max(nrow(q), nrow(nq))) %>% as.data.frame()
                        colnames(q.table) <- c("not quantified", "quantified")
                        q.table$quantified[1:nrow(q)] <- q$detectability
                        q.table$`not quantified`[1:nrow(nq)] <- nq$detectability
                  
                  save.as <- "Figure-3h.xlsx"
                  
                      wb <- createWorkbook()
                      
                      sheet <- addWorksheet(wb, "quantified")
                      writeData(wb, q, sheet=sheet, rowNames=FALSE)
                      
                      sheet <- addWorksheet(wb, "not quantified")
                      writeData(wb, nq, sheet=sheet, rowNames=FALSE)
                      
                      sheet <- addWorksheet(wb, "for Prism")
                      writeData(wb, q.table, sheet=sheet, rowNames=FALSE)
                      
                      saveWorkbook(wb, save.as, overwrite = T)
                  
                  
                  # Extended Data Figure 4g: High-detectability cysteine-containing tryptic peptides (DeepMSPeptide score > 0.5) in stereoselectively liganded proteins (protein- or cysteine-directed ABPP: log2FC in iodoacetamide-DTB labeling upon urea treatment vs. peptide quantification in cysteine-directed ABPP)
                      
                   hr.q <- sites.all %>% subset(detectability > 0.5 & enhancement.log2FC != "" & cys.abpp.detected == T, select = c(Protein, Uniprot, cys, detectability, enhancement.log2FC))
                   hr.nq <- sites.all %>% subset(detectability > 0.5 & enhancement.log2FC != "" & cys.abpp.detected == F, select = c(Protein, Uniprot, cys, detectability, enhancement.log2FC))
                  
                   hr.table <- matrix(data = "", ncol = 2, nrow = max(nrow(hr.q), nrow(hr.nq))) %>% as.data.frame()
                   colnames(hr.table) <- c("not quantified", "quantified")
                       hr.table$quantified[1:nrow(hr.q)] <- hr.q$enhancement.log2FC
                       hr.table$`not quantified`[1:nrow(hr.nq)] <- hr.nq$enhancement.log2FC
                   
                   save.as <- "Extended-Data-Figure-4g.xlsx"
                   
                       wb <- createWorkbook()
                       
                       sheet <- addWorksheet(wb, "quantified")
                       writeData(wb, hr.q, sheet=sheet, rowNames=FALSE)
                       
                       sheet <- addWorksheet(wb, "not quantified")
                       writeData(wb, hr.nq, sheet=sheet, rowNames=FALSE)
                       
                       sheet <- addWorksheet(wb, "for Prism")
                       writeData(wb, hr.table, sheet=sheet, rowNames=FALSE)
                       
                       saveWorkbook(wb, save.as, overwrite = T)

                  
              
                  
          