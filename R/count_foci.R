#' count_foci
#'
#' Count coincident foci in cropped SC and foci channel per cell
#'
#' @export count_foci
#' @param img_path, path containing image data to analyse
#' @param stage, meiosis stage of interest. Currently count_foci determines this with thresholding/ object properties in the dna channel. But will be classified using ML model in future versions.
#' @param offset_px, Pixel value offset used in thresholding of dna channel
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be small to avoid erasing foci.
#' @param foci_norm, Mean intensity to normalise all foci channels to.
#' @param annotation, Choice to output pipeline choices (recommended to knit)

#' @return foci count per cell


count_foci <- function(img_path, stage = "none", offset_px = 0.2, offset_factor = 2, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off")
{
  cell_count <- 0
  image_count <-0
  foci_counts <- 0
  antibody1_store <- 0
  antibody2_store <- 0
  if(stage == "pachytene"){
    img_path_new <- paste0(img_path,"/crops/",stage,"/")
  }
  else{
    img_path_new <- paste0(img_path,"/crops/")
  }

  file_list <- list.files(img_path_new)

  df_cols <- c("filename","cell_no","genotype","stage","foci_count", "sd_foci","mean_foci","median_foci","mean_px","median_px", "percent_on","sd_px","lone_foci")
  df_cells <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_cells) <- df_cols

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    if(stage == "pachytene"){
      filename_path_test = paste0(img_path,"/crops/",stage,"/", file)
    }
    else{
      filename_path_test = paste0(img_path,"/crops/", file)
    }

    file = filename_path_test
    if(grepl("*SYCP3.jpeg", file)){
      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(2*image, "grey")
      antibody1_store <- 1
    }
    if(grepl("*MLH3.jpeg", file)){
      file_foci = file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(  antibody1_store +antibody2_store == 2){
      antibody1_store <- 0
      antibody2_store <- 0
      cell_count <- cell_count + 1

      new_img<-img_orig
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)
      offset = offset_px
      thresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      foci_mask_crop <- img_orig_foci
      bg <- mean(img_orig_foci)
      orig_mean <- mean(img_orig_foci)
      mean_factor <- foci_norm/orig_mean
      img_orig_foci <- img_orig_foci*mean_factor

      #### normalise the foci image
      offset = offset_factor*bg
      foci_th = foci_mask_crop > bg + offset
      ### smooth it
      img_tmp_contrast = foci_mask_crop
      w = makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
      #w = makeBrush(size = 1, shape = 'gaussian', sigma = 3)
      img_flo = filter2(img_tmp_contrast, w)
      ## smooth foci channel
      foci_th = img_flo > bg + offset
      foci_label = bwlabel(foci_th)
      foci_label <- channel(foci_label, "grey")
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)

      ##### print properties of the images
      coincident_foci <- bwlabel(foci_label*strands)
      ### multiply strands by foci_label
      if(annotation == "on"){
        print("cell counter is")
        print(cell_count)
        print("original images")
        display(new_img)
        display(img_orig_foci)
        print("displaying resulting foci count")
        print("Overlay two channels")
        display(rgbImage(strands,foci_label,0*foci_label))
        print("coincident foci")
        display(colorLabels(coincident_foci))
        print("two channels, only coincident foci")
        display(rgbImage(strands,coincident_foci,coincident_foci))

      }

      overlap_no = table(coincident_foci)
      foci_per_cell <-  length(overlap_no)
      if(annotation=="on"){
        print("which counts this many foci:")
        print(foci_per_cell)

      }
      image_mat <- as.matrix(foci_mask_crop)
      image_mat <- image_mat[image_mat > 1e-06]
      #hist(image_mat)

      mean_ratio <- median(image_mat)/mean(image_mat)
      skew <- (median(image_mat)-mean(image_mat))/sd(image_mat)

      ### look at properties of the foci.
      foci_candidates <- computeFeatures.shape(foci_label)
      foci_candidates <- data.frame(foci_candidates)
      foci_areas <- foci_candidates$s.area

      ### look at properties of the overlap foci.
      overlap_candidates <- computeFeatures.shape(coincident_foci)
      overlap_candidates <- data.frame(overlap_candidates)
      overlap <- overlap_candidates$s.area

      ### number of foci NOT on an SC
      alone_foci <- nrow(foci_candidates) - foci_per_cell

      if(annotation == "on"){
        print("number of alone foci")
        print(alone_foci)
        if(alone_foci < 0){
          print("this one had a negative lone foci amount. Suspect overcounting of foci in this one:")
          display(rgbImage(strands,foci_label,0*foci_label))
        }
      }
      percent_px <- sum(overlap)/sum(foci_areas)

      if(annotation == "on"){
        print("percentage of foci channel coincident:")
        print(percent_px*100)
      }
      tryCatch({
        ### data frame stuff
        if(grepl( "++", file, fixed = TRUE) == TRUE){
          genotype <- "Fancm+/+"
        }

        if(grepl( "--", file, fixed = TRUE) == TRUE){
          genotype <- "Fancm-/-"
        }
        ### data frame stuff ends
        df_cells <- rbind(df_cells,t(c(file,cell_count,genotype,stage,foci_per_cell, sd(foci_areas),mean(foci_areas),median(foci_areas),mean(image_mat),median(image_mat),percent_px,sd(image_mat),alone_foci)))
      },
      error = function(e) {
        #what should be done in case of exception?
        str(e) # #prints structure of exception
        print("couldn't crop it")
      }
      )
    }
  }
  colnames(df_cells) <- df_cols
  return(df_cells)
}




