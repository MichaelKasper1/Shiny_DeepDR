
setwd("~/Desktop/DeepDRShiny_version2/data/")

df <- read.delim("GDSC_paper/Table_S1FG_compound_PubChem_CID_SMILES.txt")
df$PubChem_SMILES[55] <- "C1CC1COC2=CC=CC(=O)C2=C3C=C(C(=C(N3)N)C#N)C4CCNCC4"
df$filename <- df$Drug
df$filename[63] <- "VNLG_124"

library(rcdk)

for(i in 1:nrow(df)){

  if(!is.na(df$PubChem_SMILES[i])){
  smile.code.object <- parse.smiles(df$PubChem_SMILES[i])

# A temp file to save the output.
# This file will be removed later by renderImage

outfile <- paste0("GDSC_paper/DrugStructure/",df$filename[i],".png")

# Generate the png
png(outfile)
depictor <- get.depictor( width=400, height=400,zoom = 2)
img <- view.image.2d(molecule = smile.code.object[[1]],
                     depictor = depictor)

plot(NA, xlim = c(0, 2), ylim = c(0, 4), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",ann = F,bty = 'n')
rasterImage(image = img,
            xleft =  0,ybottom =  0,
            xright = 2, ytop = 4)
dev.off()
  }else{
    
    png( paste0("GDSC_paper/DrugStructure/",df$filename[i],".png"))
    par(mar = c(0, 0, 0, 0)) 
    plot(x = 0:1,                   # Create empty plot
         y = 0:1,
         ann = F,
         bty = "n",
         type = "n",
         xaxt = "n",
         yaxt = "n")
    text(x = 0.5,                   # Add text to empty plot
         y = 0.5,
         "SMILES code is not available!", 
         cex = 1.8)
    dev.off()
}



}


png( paste0("GDSC_paper/DrugStructure/Message_plot.png"))
par(mar = c(0, 0, 0, 0)) 
plot(x = 0:1,                   # Create empty plot
     y = 0:1,
     ann = F,
     bty = "n",
     type = "n",
     xaxt = "n",
     yaxt = "n")
text(x = 0.5,                   # Add text to empty plot
     y = 0.5,
     "Please select one drug from the table!", 
     cex = 1.8)
dev.off()

drug.info$filename[which(drug.info$Name == "Olaparib")] <- c("Olaparib_(rescreen)", "Olaparib")
identical(drug.info$filename, df$filename)


drug.info.out <-cbind(drug.info[,-9],df)
drug.info.out <- drug.info.out[,-9]


saveRDS(drug.info.out,"GDSC_paper/gdsc_265drug_info_update.RDS")

#Drug Info
for(i in 1:nrow(df)){
  
  if(!is.na(df$PubChem_SMILES[i])){
    smile.code.object <- parse.smiles(df$PubChem_SMILES[i])
    
    # A temp file to save the output.
    # This file will be removed later by renderImage
    
    outfile <- paste0("GDSC_paper/DrugStructure/",df$filename[i],".png")
    
    # Generate the png
    png(outfile)
    depictor <- get.depictor( width=400, height=400,zoom = 2)
    img <- view.image.2d(molecule = smile.code.object[[1]],
                         depictor = depictor)
    
    plot(NA, xlim = c(0, 2), ylim = c(0, 4), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",ann = F,bty = 'n')
    rasterImage(image = img,
                xleft =  0,ybottom =  0,
                xright = 2, ytop = 4)
    dev.off()
  }else{
    
    png( paste0("GDSC_paper/DrugStructure/",df$filename[i],".png"))
    par(mar = c(0, 0, 0, 0)) 
    plot(x = 0:1,                   # Create empty plot
         y = 0:1,
         ann = F,
         bty = "n",
         type = "n",
         xaxt = "n",
         yaxt = "n")
    text(x = 0.5,                   # Add text to empty plot
         y = 0.5,
         "SMILES code is not available!", 
         cex = 1.8)
    dev.off()
  }
  
  
  
}
