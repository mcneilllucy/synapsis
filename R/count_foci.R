#' count_foci
#'
#' Calculates coincident foci in synaptonemal complex and foci channel, per cell
#'
#' In this function, masks for the synaptonemal complex (SC) and foci channel
#' are created from the saved crops of single/individual cells.
#' These masks are computed using (optional) input parameters related
#' to meiosis stage/ how well spread chromosomes are (for the former)
#' and related to smoothing, thresholding and how "crowded" foci are for the
#' latter. Finally, these two masks are multiplied, and the number of
#' objects found with EBImage's computeFeatures are the colocalizing foci.
#'
#' The file, cell number, foci count etc. are output as a data frame.
#'
#' @export count_foci
#' @param img_path, path containing crops to analyse
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the synaptonemal complex channel
#' by previosly calling the get_pachytene function.
#' Note that if using this option, the count_foci function requires that the
#' input directory contains a folder called “pachytene” with the crops in it.
#' @param offset_px, Pixel value offset used in thresholding of synaptonemal complex channel
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small
#' to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @param foci_norm, Mean intensity to normalise all foci channels to.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param channel1_string String appended to the files showing the channel
#' illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel
#' illuminating synaptonemal complexes. Defaults to SYCP3
#' @param file_ext file extension of your images e.g. tiff jpeg or png.
#' @param KO_str string in filename corresponding to knockout genotype.
#' Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype.
#' Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout.
#' Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout.
#' Defaults to +/+.
#' @param watershed_stop Stop default watershed method with "on"
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @param crowded_foci TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have foci > 100 or so.
#' @param artificial_amp_factor Amplification of foci channel, for annotation only.
#' @param strand_amp multiplication of strand channel to make masks
#' @param min_foci minimum pixel area for a foci. Depends on your dpi etc. Defaults to 4
#' @param disc_size size of disc for local background calculation in synaptonemal complex channel
#' @param modify_problematic option for synapsis to try and "save" images which
#' have likely been counted incorrectly due to a number of reasons. Default
#' settings are optimized for mouse pachytene. Defaults to "off"
#' @param disc_size_foci size of disc for local background calculation in foci channel
#' @param C_weigh_foci_number choose crispness criteria- defaults to TRUE to use
#' C1 (weighing with number). Otherwise set to FALSE to use C2
#' @param C1 Default crispness criteria = sd(foci_area)/(mean(foci_area)+1)
#' @param C2 Alternative crisp criteria.
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' foci_counts <- count_foci(demo_path,offset_factor = 3, brush_size = 3,
#' brush_sigma = 3, annotation = "on",stage = "pachytene")
#' @author Lucy McNeill
#' @return data frame with foci count per cell


count_foci <- function(img_path, stage = "none", offset_px = 0.2, offset_factor = 2, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off",channel2_string = "SYCP3", channel1_string = "MLH3",file_ext = "jpeg", KO_str = "--",WT_str = "++",KO_out = "-/-", WT_out = "+/+", watershed_stop = "off", watershed_radius = 1, watershed_tol = 0.05, crowded_foci = TRUE, artificial_amp_factor = 1, strand_amp = 2, min_foci =-1, disc_size = 51, modify_problematic = "off", disc_size_foci = 5, C1 = 0.02, C2 = 0.46, C_weigh_foci_number = TRUE)
{
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  if(stage == "pachytene"){
    img_path_new <- paste0(img_path,"/crops/",stage,"/")
  }
  else{
    img_path_new <- paste0(img_path,"/crops/")
  }
  file_list <- list.files(img_path_new)
  df_cols <- c("filename","cell_no","genotype","stage","foci_count", "sd_foci","mean_foci","median_foci","mean_px","median_px", "percent_on","sd_px","lone_foci","comments","verdict","C1")
  df_cells <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_cells) <- df_cols
  for (img_file in file_list){
    if(stage == "pachytene"){
      filename_path_test <- paste0(img_path,"/crops/",stage,"/", img_file)
    }
    else{
      filename_path_test <- paste0(img_path,"/crops/", img_file)
    }
    img_file <- filename_path_test
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
      file_sc <- img_file
      image_count <- image_count +1
      image <- readImage(file_sc)
      img_orig <- channel(strand_amp*image, "grey")
      antibody1_store <- 1
    }
    if(grepl(paste0('*',channel1_string,'.',file_ext,'$'), img_file)){
      file_foci <- img_file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      antibody2_store <- 1
    }
    if(antibody1_store +antibody2_store == 2){
      C1_search <- TRUE
      antibody1_store <- 0
      antibody2_store <- 0
      discrepant_category <- 0
      cell_count <- cell_count + 1
      df_cells <- get_coincident_foci(offset_px, offset_factor, brush_size,
                                      brush_sigma, annotation,watershed_stop,
                                      watershed_radius, watershed_tol,
                                      crowded_foci, artificial_amp_factor,
                                      strand_amp, disc_size, disc_size_foci,
                                      img_file,cell_count,img_orig,
                                      img_orig_foci, stage,
                                      WT_str,KO_str,WT_out,KO_out, C1_search,
                                      discrepant_category,C1,C2,df_cells, C_weigh_foci_number)
    }
  }
  colnames(df_cells) <- df_cols
  return(df_cells)
}



#' annotate_foci_counting
#'
#' Contains all plotting routines for count foci annotation
#'
#' @param img_file cell's file name
#' @param cell_count unique cell counter
#' @param img_orig original strand crop
#' @param img_orig_foci cropped foci channel
#' @param artificial_amp_factor amplification factor
#' @param strands black white mask of strand channel
#' @param coincident_foci mask of overlap between strand and foci channel
#' @param foci_label black and white mask of foci channel
#' @param percent_px percentage of foci mask that coincides with strand channel
#' small number indicates potentially problematic image.
#' @param alone_foci estimated number of foci that are NOT on a strand.
#' @param foci_per_cell number of foci counted per cell
#' @return displays key steps from raw image to coincident foci count
#'
annotate_foci_counting <- function(img_file,cell_count,img_orig,img_orig_foci,
                                   artificial_amp_factor,strands,
                                   coincident_foci, foci_label,alone_foci,
                                   percent_px,foci_per_cell){
  cat("\n at file",img_file, sep = " ")
  cat("\n cell counter is", cell_count, sep= " ")
  cat("\n original images")
  plot(img_orig)
  plot(img_orig_foci)
  ch1 <-channel(img_orig,"grey")
  ch2 <- channel(artificial_amp_factor*img_orig_foci,"grey")
  bluered <- rgbImage(ch1, ch2, 0*ch1)
  cat("\n displaying resulting foci count plots. Overlay two channels:")
  plot(rgbImage(ch1,ch2,0*img_orig))
  cat("\n displaying resulting masks. Overlay two masks:")
  plot(rgbImage(strands,foci_label,0*img_orig))
  cat("\n coincident foci:")
  plot(colorLabels(coincident_foci))
  cat("\n two channels, only coincident foci")
  plot(rgbImage(strands,coincident_foci,coincident_foci))
  cat("\n which counts this many foci:",foci_per_cell, sep = " ")
  cat("\n number of alone foci",alone_foci, sep = " ")
  if(alone_foci < 0){
    cat("\n this one had a negative lone foci amount. Suspect overcounting of foci in this one:")
    plot(rgbImage(strands,foci_label,0*foci_label))
    #alone_foci <- 0
  }
  tryCatch({
    if(percent_px < 0.3){
      cat("\n file name", img_file, "coincides", percent_px*100,
          "percent. So decided to reduce the local background thresholding
              requirements.", sep = " ")
    }
  },
  error = function(e) {
    #what should be done in case of exception?
    cat("\n percent on was a NaN. Suspect no/low signal")
  }
  )
}

#' append_data_frame
#'
#' applies new row to data frame
#'
#' @param img_file cell's file name
#' @param cell_count unique cell counter
#' @param KO_str string in filename corresponding to knockout genotype.
#' Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype.
#' Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout.
#' Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout.
#' Defaults to +/+.
#' @param foci_areas pixel area of each foci
#' @param df_cells current data frame
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the synaptonemal complex channel
#' by previosly calling the get_pachytene function.
#' Note that if using this option, the count_foci function requires that the
#' input directory contains a folder called “pachytene” with the crops in it.
#' @param foci_per_cell foci count for cell
#' @param image_mat matrix with all pixel values above zero
#' @param percent_px percentage of foci mask that coincides with strand channel
#' small number indicates potentially problematic image.
#' @param alone_foci estimated number of foci that are NOT on a strand.
#' @param discrepant_category estimated number of foci that are NOT on a strand.
#' @param C1 criteria
#'
#' @return data frame with new row
#'
append_data_frame <- function(WT_str,KO_str,WT_out,KO_out,img_file,foci_areas,df_cells,cell_count,stage,foci_per_cell,image_mat,percent_px,alone_foci,discrepant_category, C1){
  tryCatch({
    ### data frame stuff
    if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
      genotype <- WT_out
    }
    if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
      genotype <- KO_out
    }

    if(foci_per_cell < 2){
      ### force statistics about foci areas to all be zero rather than fail
      foci_areas <- c(0,0)
    }
    ### data frame stuff ends
    if(discrepant_category == 0){
      discrepant_comment <- "meets crisp-ness criteria"
      verdict <- "keep"

    }
    else if(discrepant_category ==1){
      discrepant_comment <- "does not meet crisp-ness criteria"
      verdict <- "discard"
    }
    else if(discrepant_category >1){
      discrepant_comment <- "does not meet crisp-ness criteria"
      verdict <- "discard"
    }

    df_cells <- rbind(df_cells,t(c(img_file,cell_count,genotype,stage,
                                   foci_per_cell, sd(foci_areas),mean(foci_areas),
                                   median(foci_areas),mean(image_mat),median(image_mat),
                                   percent_px,sd(image_mat),alone_foci,discrepant_comment,
                                   verdict,C1)))
  },
  error = function(e) {
    #what should be done in case of exception?
    str(e) # #prints structure of exception
    cat("\n something went wrong while making the data frame")
  }
  )
  return(df_cells)
}


### annotate when the flag has been raised that the percent of foci on is too small.
## Nan case indicates zero, so should move to a separate case of dim foci possibly.
## Change the thresholding to longer be local background (but still smoothed)
## perhaps base the new threshold value on the mean pixel value... rather then difference compared to background?
### then make a test about whether there was an improvement....


#### annonate when a flag has been raised about the XY stuff:
### if sd(foci) is >10* typical value, then you have an XY blob there. Remove it from the mask?

#' remove_XY
#'
#' applies new row to data frame
#'
#' @param foci_label black and white mask of foci channel
#' @param foci_candidates computeFeatures data frame of foci channel
#' @param foci_areas the areas of the foci objects
#'
#'
#'
#' @return mask with XY blob removed
#'
#'
remove_XY <- function(foci_label, foci_candidates, foci_areas){
  ### if the standard deviation was too big
  ### loop over and delete the giant blobs.
  ### after this function is run, coincident foci needs to be computed again.
  median_px <- median(foci_areas)
  #### takes in foci mask. Outputs the new mask without giant blob.
  x <- computeFeatures.shape(foci_label)
  x <- data.frame(x)
  OOI <- nrow(x)
  counter <- 0
  retained <- foci_label
  while(counter<OOI){
    counter <- counter+1
    pixel_area <- x$s.area[counter]
    # if statement checking if it's the big blob
    if(pixel_area> 10*median_px){
      retained <- as.numeric(retained)*rmObjects(foci_label, counter, reenumerate = TRUE)
    }
  }
  return(retained)
  ## multiply with the original foci_label

}

#' make_foci_mask
#'
#' creates foci mask for foci channel crop
#'
#' @param offset_factor Pixel value offset used in thresholding of foci channel
#' @param bg background value- currently just mean pixel value of whole image
#' @param crowded_foci TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have foci > 100 or so.
#' @param img_orig_foci cropped foci channel
#' @param brush_size size of brush to smooth the foci channel. Should be small
#' to avoid erasing foci.
#' @param brush_sigma sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @param disc_size_foci size of disc for local background calculation in foci channel
#' @return foci mask
#'
#'
make_foci_mask <- function(offset_factor,bg,crowded_foci,img_orig_foci,
                           brush_size,brush_sigma,disc_size_foci){
  foci_mask_crop <- img_orig_foci
  offset <- offset_factor*bg
  if(crowded_foci == TRUE){
    ### do local background calculation instead. but still don't smooth.
    #foci_th <- foci_mask_crop > bg + offset
    ### using local bg
    disc_size <- disc_size_foci
    new_img<-foci_mask_crop
    disc <- makeBrush(disc_size, "disc")
    disc <- disc / sum(disc)
    localBackground <- filter2(new_img, disc)
    offset <- 1/(2.8*offset_factor)
    foci_th <- (new_img - localBackground > offset)
  }
  else{
    ### smooth it
    img_tmp_contrast <- foci_mask_crop
    w <- makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
    img_flo <- filter2(img_tmp_contrast, w)
    ### using local bg
    disc_size <- disc_size_foci
    new_img<-img_flo
    disc <- makeBrush(disc_size, "disc")
    disc <- disc / sum(disc)
    localBackground <- filter2(new_img, disc)
    offset <- 1/(2.8*offset_factor)
    foci_th <- (new_img - localBackground > offset)

  }
  foci_label <- bwlabel(foci_th)
  return(foci_label)
}

#' make_strand_mask
#'
#' creates strand mask for strand channel crop
#'
#' @param offset_px, Pixel value offset used in thresholding of synaptonemal complex channel
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the synaptonemal complex channel
#' by previosly calling the get_pachytene function.
#' Note that if using this option, the count_foci function requires that the
#' input directory contains a folder called “pachytene” with the crops in it.
#' @param img_orig original strand crop
#' @param disc_size size of disc for local background calculation in synaptonemal complex channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small
#' to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @return strand mask
#'
make_strand_mask <- function(offset_px, stage, img_orig, disc_size,brush_size,brush_sigma){
  img_tmp_contrast <- img_orig
  w <- makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
  img_flo <- filter2(img_orig, w)
  new_img<-img_flo
  disc <- makeBrush(disc_size, "disc")
  disc <- disc / sum(disc)
  localBackground <- filter2(new_img, disc)
  offset <- offset_px
  if(stage == "pachytene"){
    thresh_crop <- (new_img - localBackground > offset)

  }
  else{
    thresh_crop <- new_img > offset
  }
  strands <- bwlabel(thresh_crop)
  return(strands)
}

#' get_overlap_mask
#'
#' creates mask for coincident foci
#'
#' @param strands black white mask of strand channel
#' @param watershed_stop Stop default watershed method with "on"
#' @param foci_label black and white mask of foci channel
#' @param img_orig_foci cropped foci channel
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @return mask with coincident foci on strands
#'
get_overlap_mask<- function(strands, foci_label, watershed_stop, img_orig_foci, watershed_radius, watershed_tol){
  num_strands <- computeFeatures.shape(strands)
  num_strands <- data.frame(num_strands)
  if(watershed_stop != "off"){
    coincident_foci <- bwlabel(foci_label*strands)
  }
  else{
    foci_th <- channel(foci_label,"grey")
    coincident_foci <- watershed(bwlabel(foci_th*strands)*as.matrix(img_orig_foci),tolerance=watershed_tol, ext=watershed_radius)
  }
  return(coincident_foci)
}

#' get_foci_per_cell
#'
#' creates mask for coincident foci
#'
#' @param img_file cell's file name
#' @param offset_px, Pixel value offset used in thresholding of synaptonemal complex channel
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the synaptonemal complex channel
#' by previosly calling the get_pachytene function.
#' Note that if using this option, the count_foci function requires that the
#' input directory contains a folder called “pachytene” with the crops in it.
#' @param strands black white mask of strand channel
#' @param watershed_stop Stop default watershed method with "on"
#' @param foci_label black and white mask of foci channel
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param cell_count unique cell counter
#' @param img_orig original strand crop
#' @param img_orig_foci cropped foci channel
#' @param artificial_amp_factor amplification factor
#' @param coincident_foci mask of coincident foci
#'
#' @return number of foci per cell
#'
get_foci_per_cell <- function(img_file,offset_px,stage,strands,watershed_stop,foci_label, annotation, cell_count, img_orig, img_orig_foci, artificial_amp_factor, coincident_foci){
  coincident_df <- data.frame(computeFeatures.shape(coincident_foci))
  if(annotation == "on"){
  }
  coincident_df <- coincident_df[coincident_df$s.area,]
  ### multiply strands by foci_label
  foci_per_cell <- nrow(coincident_df)
  return(foci_per_cell)
}


#' annotate_foci_counting_adjusted
#'
#' Contains all plotting routines for count foci annotation
#'
#' @param img_file cell's file name
#' @param cell_count unique cell counter
#' @param img_orig original strand crop
#' @param img_orig_foci cropped foci channel
#' @param artificial_amp_factor amplification factor
#' @param strands black white mask of strand channel
#' @param coincident_foci mask of overlap between strand and foci channel
#' @param foci_label black and white mask of foci channel
#' @param percent_px percentage of foci mask that coincides with strand channel
#' small number indicates potentially problematic image.
#' @param alone_foci estimated number of foci that are NOT on a strand.
#' @param foci_per_cell number of foci counted per cell
#' @return displays key steps from raw image to coincident foci count
#'
annotate_foci_counting_adjusted <- function(img_file,cell_count,img_orig,img_orig_foci,artificial_amp_factor,strands,coincident_foci, foci_label,alone_foci,percent_px,foci_per_cell){
  cat("\n at file",img_file, sep = " ")
  cat("\n cell counter is", cell_count, sep= " ")
  cat("\n original images")
  plot(img_orig)
  plot(img_orig_foci)
  ch1 <-channel(img_orig,"grey")
  ch2 <- channel(artificial_amp_factor*img_orig_foci,"grey")
  bluered <- rgbImage(ch1, ch2, 0*ch1)
  cat("\n displaying resulting foci count plots. Overlay two channels:")
  plot(rgbImage(ch1,ch2,0*img_orig))
  cat("\n displaying resulting masks. Overlay two masks:")
  plot(rgbImage(strands,foci_label,0*img_orig))
  cat("\n coincident foci:")
  plot(colorLabels(coincident_foci))
  cat("\n two channels, only coincident foci")
  plot(rgbImage(strands,coincident_foci,coincident_foci))
  cat("\n which counts this many foci:",foci_per_cell, sep = " ")
  cat("\n number of alone foci",alone_foci, sep = " ")
}


#' get_C1
#'
#' calculates the statistic to compare to crisp_criteria, which determines
#' whether the foci count will be reliable
#'
#' @param C_weigh_foci_number choose crispness criteria- defaults to TRUE to use
#' C1 (weighing with number). Otherwise set to FALSE to use C2
#' @param foci_areas pixel area of each foci
#' @param foci_per_cell foci count for cell
#' @return statistic to comapre to crisp_criteria

get_C1 <- function(foci_areas, foci_per_cell, C_weigh_foci_number){
  if(C_weigh_foci_number == TRUE){
    C1 <- sd(foci_areas)/((mean(foci_areas))*(foci_per_cell+1))
  }
  else{
    C1 <- sd(foci_areas)/((mean(foci_areas))*(foci_per_cell+1)/(foci_per_cell+1))
  }
  return(C1)
}


#' get_coincident_foci
#'
#' calculates the statistic to compare to crisp_criteria, which determines
#' whether the foci count will be reliable
#' @param img_file cell's file name
#' @param cell_count unique cell counter
#' @param img_orig original strand crop
#' @param img_orig_foci cropped foci channel
#' @param C1_search TRUE or FALSE whether the image is still being modified
#' until it meets the crispness criteria
#' @param discrepant_category estimated number of foci that are NOT on a strand.
#' @param df_cells current data frame
#' @param C_weigh_foci_number choose crispness criteria- defaults to TRUE to use
#' C1 (weighing with number). Otherwise set to FALSE to use C2
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the synaptonemal complex channel
#' by previosly calling the get_pachytene function.
#' Note that if using this option, the count_foci function requires that the
#' input directory contains a folder called “pachytene” with the crops in it.
#' @param offset_px, Pixel value offset used in thresholding of synaptonemal complex channel
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small
#' to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param KO_str string in filename corresponding to knockout genotype.
#' Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype.
#' Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout.
#' Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout.
#' Defaults to +/+.
#' @param watershed_stop Stop default watershed method with "on"
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @param crowded_foci TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have foci > 100 or so.
#' @param artificial_amp_factor Amplification of foci channel, for annotation only.
#' @param strand_amp multiplication of strand channel to make masks
#' @param disc_size size of disc for local background calculation in synaptonemal complex channel
#' @param disc_size_foci size of disc for local background calculation in foci channel
#' @param C_weigh_foci_number choose crispness criteria- defaults to TRUE to use
#' C1 (weighing with number). Otherwise set to FALSE to use C2
#' @param C1 Default crispness criteria = sd(foci_area)/(mean(foci_area)+1)
#' @param C2 Alternative crisp criteria.
#'
#' @return data frame with new row with most recent foci per cell appended
#'
#'
get_coincident_foci <- function(offset_px, offset_factor, brush_size, brush_sigma, annotation, watershed_stop, watershed_radius, watershed_tol, crowded_foci, artificial_amp_factor, strand_amp, disc_size, disc_size_foci,img_file,cell_count,img_orig,img_orig_foci,stage,WT_str,KO_str,WT_out,KO_out, C1_search,discrepant_category,C1,C2,df_cells,C_weigh_foci_number){
  strands <- make_strand_mask(offset_px, stage, img_orig, disc_size,brush_size,brush_sigma)
  color_img_strands<- colorLabels(strands, normalize = TRUE)
  num_strands <- computeFeatures.shape(strands)
  num_strands <- data.frame(num_strands)
  bg <- mean(img_orig_foci)
  if(C_weigh_foci_number == TRUE){
    crisp_criteria <- C1
  }
  else{
    crisp_criteria <- C2
  }
  while(C1_search == TRUE && discrepant_category < 2){
    foci_label <- make_foci_mask(offset_factor,bg,crowded_foci,img_orig_foci,brush_size,brush_sigma,disc_size_foci)
    foci_label <- channel(foci_label, "grey")
    foci_candidates <- data.frame(computeFeatures.shape(foci_label))
    foci_areas <- foci_candidates$s.area
    coincident_foci <- get_overlap_mask(strands, foci_label, watershed_stop, img_orig_foci, watershed_radius, watershed_tol)
    overlap_candidates <- data.frame(computeFeatures.shape(coincident_foci))
    overlap <- overlap_candidates$s.area
    percent_px <- sum(overlap)/sum(foci_areas)
    foci_per_cell <- get_foci_per_cell(img_file,offset_px,stage,strands,watershed_stop,foci_label,annotation, cell_count,img_orig, img_orig_foci, artificial_amp_factor,coincident_foci)
    C1 <- get_C1(foci_areas, foci_per_cell,C_weigh_foci_number)
    if(is.na(C1)){
      C1 <- 100
    }
    if(C1 < crisp_criteria){
      C1_search <- FALSE
    }
    else{
      discrepant_category <- discrepant_category + 1
      C1_search <- FALSE
    }
  }
  alone_foci <- nrow(foci_candidates) - foci_per_cell
  if(annotation=="on"){
    annotate_foci_counting(img_file,cell_count,img_orig,img_orig_foci,artificial_amp_factor,strands,coincident_foci, foci_label,alone_foci,percent_px,foci_per_cell)
  }
  image_mat <- as.matrix(img_orig_foci)
  image_mat <- image_mat[image_mat > 1e-06]
  df_cells <- append_data_frame(WT_str,KO_str,WT_out,KO_out,img_file,foci_areas,df_cells,cell_count,stage,foci_per_cell,image_mat,percent_px,alone_foci,discrepant_category,C1)
  return(df_cells)
  }
