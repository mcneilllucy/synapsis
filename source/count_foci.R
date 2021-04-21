count_foci <- function(file_list, img_path)
{
  # input :
  
  # output : a bunch of output jpegs? Or save them all?
  
  
  
  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  pair <- 0
  foci_counts <- 0
  
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
      
      #print("I have a pair of images")
      #display(img_orig)
      #display(img_orig_foci)
      
      new_img<-img_orig
      #### now see which have the right amount of strands 
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)
      offset = 0.1
      nucBadThresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(nucBadThresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      ### define a foci mask
      w = makeBrush(size = 51, shape = 'gaussian', sigma = 10)
      img_flo_foci = filter2(img_orig_foci, w)
      ## do adaptive thresholding on the image.. this should solve the issue...
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(img_orig_foci, disc)
      offset = 0.1
      nucBadThresh_foci = (img_orig_foci - localBackground > offset)
      foci_mask <- bwlabel(nucBadThresh_foci)
      #display(foci_mask)
      
      foci_mask_crop <- foci_mask
      
      
      bg <- mean(img_orig_foci)
      if(bg < 0.05 && bg > 0.02){
        display(foci_mask)
        offset = 4*bg
        foci_th = foci_mask_crop > bg + offset
        ### smooth it
        ### maybe up the contrast first??
        img_tmp_contrast = foci_mask_crop
        print("cell counter is")
        print(cell_count)
        display(foci_mask_crop)
        w = makeBrush(size = 3, shape = 'gaussian', sigma = 1)
        img_flo = filter2(img_tmp_contrast, w)
        ## only choose objects above bright pixel value
        #foci_th = foci_mask_crop > 0.2
        #foci_th = img_flo > 0.2
        
        #display(foci_th)
        foci_label = bwlabel(foci_th)
        foci_label <- channel(foci_label, "grey")
        display(strands)
        display(foci_label)
        
        ### multiply strands by foci_label
        display(rgbImage(strands,foci_label,foci_label))
        coincident_foci <- bwlabel(foci_label*strands)
        display(colorLabels(coincident_foci))
        overlap_no = table(coincident_foci)
        foci_per_cell <-  length(overlap_no)
        #print(file)
        #print("foci count was:")
        #print(foci_per_cell)
        if (foci_per_cell > 0){
          foci_counts <- append(foci_counts,foci_per_cell)
        }
        
        
        
      }
      
      
      ###
    }
  
  }
  hist(foci_counts)
  print(mean(foci_counts))
  print(median(foci_counts))
  print(sd(foci_counts))
  return(foci_counts)

}
