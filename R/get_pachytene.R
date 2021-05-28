#' get_pachytene
#'
#' Identifies crops in pachytene
#'
#' @import EBImage
#' @export
#' @param file_list The file list
#' @param img_path The path
#' @return Pairs of foci and SC channel crops for pachytene


get_pachytene <- function(img_path)
{
  # input :

  # output : a bunch of output jpegs? Or save them all?



  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  pachytene_count <- 0
  setwd(img_path)

  img_path_new <- paste0(img_path,"/crops/")
  print(img_path_new)

  setwd(img_path_new)
  dir.create("pachytene")
  file_list <- list.files(img_path_new)
  print(file_list)

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    setwd(img_path_new)
    if(grepl("*SYCP3.jpeg", file)){
      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(image, "grey")
      antibody1_store <- 1
    }
    if(grepl("*MLH3.jpeg", file)){
      file_foci = file
      #print(file_foci)
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(antibody1_store + antibody2_store ==2){
      antibody1_store <- 0
      antibody2_store <- 0
      print("I have a pair")
      new_img<-img_orig
      #### now see which have the right amount of strands
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)

      # lax
       offset = 0.1

      # strict
      #offset = 0.2
      thresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      print(nrow(num_strands))
      #### segment the strands
      if (nrow(num_strands)<22 && nrow(num_strands)>5){
        ### identified a good image. count foci
        display(new_img)
        display(strands)
        pachytene_count <- pachytene_count + 1

        file_dna <- tools::file_path_sans_ext(file_dna)
        filename_crop = paste0("./pachytene/", file_dna,".jpeg")
        writeImage(img_orig, filename_crop)

        file_foci <- tools::file_path_sans_ext(file_foci)
        filename_crop_foci = paste0("./pachytene/", file_foci, ".jpeg")
        writeImage(img_orig_foci, filename_crop_foci)

        print("two filenames are")
        print(filename_crop)
        print(filename_crop_foci)
        print("end filename")

      }
      ###
    }
  }
print("number of cells kept")
print(pachytene_count)
}

