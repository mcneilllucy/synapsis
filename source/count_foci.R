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
      img_orig <- channel(2*image, "grey")
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
      display(new_img)
      #### now see which have the right amount of strands
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)
      offset = 0.2
      nucBadThresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(nucBadThresh_crop)
      display(strands)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      foci_mask_crop <- img_orig_foci
      bg <- mean(img_orig_foci)
      orig_mean <- mean(img_orig_foci)
      mean_factor <- 0.01/orig_mean
      img_orig_foci <- img_orig_foci*mean_factor
      display(img_orig_foci)
      #### normalise the foci image
      if(bg < 1 && bg > 0){
        offset = 2*bg
        foci_th = foci_mask_crop > bg + offset
        ### smooth it
        ### maybe up the contrast first??
        img_tmp_contrast = foci_mask_crop

        print("cell counter is")
        print(cell_count)
        #display(foci_mask_crop)
        w = makeBrush(size = 3, shape = 'gaussian', sigma = 3)
        img_flo = filter2(img_tmp_contrast, w)
        ## only choose objects above bright pixel value
        #foci_th = foci_mask_crop > 0.2
        #foci_th = img_flo > 0.2
        display(img_flo)

        ## smooth foci channel
        foci_th = img_flo > bg + offset
        #foci_th = img_flo > 0.05

        display(foci_th)
        foci_label = bwlabel(foci_th)
        foci_label <- channel(foci_label, "grey")
        display(colorLabels(strands))
        num_strands <- computeFeatures.shape(strands)
        num_strands <- data.frame(num_strands)

        ##### print properties of the images

        ### multiply strands by foci_label
        display(rgbImage(strands,foci_label,0*foci_label))
        coincident_foci <- bwlabel(foci_label*strands)
        display(colorLabels(coincident_foci))
        overlap_no = table(coincident_foci)
        foci_per_cell <-  length(overlap_no)
        print(foci_per_cell)
        #print(file)

        image_mat <- as.matrix(foci_mask_crop)
        image_mat <- image_mat[image_mat > 1e-06]
        hist(image_mat)

        mean_ratio <- median(image_mat)/mean(image_mat)
        skew <- (median(image_mat)-mean(image_mat))/sd(image_mat)

        ### look at properties of the foci.
        foci_candidates <- computeFeatures.shape(foci_label)
        foci_candidates <- data.frame(foci_candidates)
        foci_areas <- foci_candidates$s.area

        if (sd(foci_areas)<20 && foci_per_cell >0){
          foci_counts <- append(foci_counts,foci_per_cell)
        }

      }


      ###
    }

  }
  hist(foci_counts, breaks =7 )
  print(mean(foci_counts))
  print(median(foci_counts))
  print(sd(foci_counts))
  return(foci_counts)

}
