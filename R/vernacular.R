# Get norwegian and english common vernacular names
# English names are taken from the gbif backbone with the rgbif package
# Norwegian names are taken from a list of species names from Artsdatabanken. 
# The translations are written to file. 
# This should be done after running the fitModels script. 
# The file is included in the shiny bundle and loaded on initiation.
# That way one doesn't need to wait for the name_lookup functin on initiation as that takes 15-20 seconds,
# and may not give the same result each time.
# There are several other possible rutes to take. Taxize is one that I didn't have great success with.


library(rgbif)
library(dplyr)
library(stringr)


# first a test


ver <- name_lookup(query='Ajuga reptans', rank="species", return="names")
ver <- dplyr::bind_rows(ver, .id = "column_label")
ver2 <- ver$vernacularName[ver$language == "eng"]
ver2 <- ver2[!is.na(ver2)]
 # 26 common names found - which one to chose?
 # just get the first one
ver2 <- ver2[1] # Common bugle


# Now lets make it automatic for a lenger list if species
myS <- list.files("sdmModels/", pattern = ".sdm")
myS2 <- as.list(NA)
for(i in 1:length(myS)){
  myS2[i] <- 
    paste(
      stringr::str_split(myS[i], "_")[[1]][1],
      stringr::str_split(myS[i], "_")[[1]][2],
      collapse = " ")
}
myS3x <- unique(as.character(myS2))

# Get norwegian names
namelist <- readr::read_delim("Artsnavnebase_fork.csv", 
                              "\t", escape_double = FALSE, locale = readr::locale(encoding = "ISO-8859-1"),
                              trim_ws = TRUE)
namelist$sp <- paste(namelist$Slekt, namelist$Art)
namelist <- namelist[!duplicated(namelist$sp),]  # 13.5k names
myS3 <- data.frame(myS3 = myS3x)
for(i in 1:nrow(myS3)){
  myS3[i,2] <-  namelist$PopulærnavnBokmål[grep(myS3x[i], namelist$sp, fixed = TRUE)]
}
# Manually fixing erroneous translations (see issue #8)
myS3$V2[myS3$myS3=="Arnica montana"]      <- "solblom"
myS3$V2[myS3$myS3=="Carex lepidocarpa"]   <- "nebbstarr"
myS3$V2[myS3$myS3=="Lathyrus palustris"]  <- "myrflatbelg"
myS3$V2[myS3$myS3=="Malus sylvestris"]    <- "villeple"
myS3$V2[myS3$myS3=="Pseudorchis albida"]  <- "hvitkurle"
myS3$V2[myS3$myS3=="Thalictrum simplex"]  <- "rankfrøstjerne"
myS3$V2[myS3$myS3=="Thymus praecox"]      <- "kryptimian"


sci <- as.character(myS3$myS3)
vern_nor <- as.character(myS3$V2)
rm(myS2, namelist, myS, myS3x, myS3)

vern_eng <- NA
for (i in 1:length(sci)){
  s <- sci[i]
  ver <- name_lookup(query=s, rank="species", return="names")
  ver <- dplyr::bind_rows(ver, .id = "column_label")
  ver2 <- ver$vernacularName[ver$language == "eng"]
  ver2 <- ver2[!is.na(ver2)]
  if(length(ver2)>0){
    vern_eng[i] <- tolower(ver2)
  } else{
    vern_eng[i] <- NA
  }
}

rm(i, s, ver, ver2)


vern_eng[is.na(vern_eng)] <- sci[is.na(vern_eng)]

#View(cbind(vern_eng, vern_eng2, sci))
vern_eng2[is.na(vern_eng2)] <- "English name missing"

namelist <- data.frame(sci, vern_nor, vern_eng, vern_eng2)

saveRDS(namelist, 'namelist.RData')
