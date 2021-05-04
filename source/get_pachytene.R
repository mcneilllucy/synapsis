get_pachytene <- function(file_list, img_path)
{
  # input :

  # output : a bunch of output jpegs? Or save them all?



  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  pair <- 0
  pachytene_count <- 0

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    setwd(img_path)
    if(grepl("*dna.jpeg$", file)){
      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(image, "grey")
      pair <- 0
    }
    if(grepl("*foci.jpeg$", file)){
      file_foci = file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      pair <- 1
    }
    if(pair ==1){
      new_img<-img_orig
      #### now see which have the right amount of strands
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)

      # lax
       offset = 0.1

      # strict
      #offset = 0.2
      nucBadThresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(nucBadThresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      #### segment the strands
      if (width(num_strands)<22 && width(num_strands)>5){
        ### identified a good image. count foci
        display(new_img)
        display(strands)
        pachytene_count <- pachytene_count + 1

        filename_crop = paste0("./pachytene/", file,"-dna.jpeg")
        writeImage(img_orig, filename_crop)

        filename_crop_foci = paste0("./pachytene/", file, "-foci.jpeg")
        writeImage(img_orig_foci, filename_crop_foci)

      }


      ###
    }

  }
print("number of cells kept")
print(pachytene_count)
}
