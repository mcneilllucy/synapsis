#' count_foci
#'
#' Calculates coincident foci in SC and foci channel, per cell
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
#' this with thresholding/ object properties in the dna channel. But will be
#' classified using ML model in future versions.
#' @param offset_px, Pixel value offset used in thresholding of dna channel
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
#' @param watershed_stop Turn off default watershed method with "off"
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @param crowded_foci TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have foci > 100 or so.
#' @param artificial_amp_factor Amplification of foci channel, for annotation only.
#' @param strand_amp multiplication of strand channel to make masks
#' @param min_foci minimum pixel area for a foci. Depends on your dpi etc. Defaults to 4
#' @param disc_size size of disc for local background calculation in dna channel
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' foci_counts <- count_foci(demo_path,offset_factor = 3, brush_size = 3,
#' brush_sigma = 3, annotation = "on",stage = "pachytene")
#' @author Lucy McNeill
#' @return foci count per cell


count_foci <- function(img_path, stage = "none", offset_px = 0.2, offset_factor = 2, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off",channel2_string = "SYCP3", channel1_string = "MLH3",file_ext = "jpeg", KO_str = "--",WT_str = "++",KO_out = "-/-", WT_out = "+/+", watershed_stop = "off", watershed_radius = 1, watershed_tol = 0.05, crowded_foci = TRUE, artificial_amp_factor = 1, strand_amp = 2, min_foci =-1, disc_size = 51)
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
  df_cols <- c("filename","cell_no","genotype","stage","foci_count", "sd_foci","mean_foci","median_foci","mean_px","median_px", "percent_on","sd_px","lone_foci")
  df_cells <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_cells) <- df_cols
  ## for each image that is *-dna.jpeg,
  for (img_file in file_list){
    if(stage == "pachytene"){
      filename_path_test <- paste0(img_path,"/crops/",stage,"/", img_file)
    }
    else{
      filename_path_test <- paste0(img_path,"/crops/", img_file)
    }
    img_file <- filename_path_test
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
      file_dna <- img_file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(strand_amp*image, "grey")
      antibody1_store <- 1
    }
    if(grepl(paste0('*',channel1_string,'.',file_ext,'$'), img_file)){
      file_foci <- img_file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(antibody1_store +antibody2_store == 2){
      antibody1_store <- 0
      antibody2_store <- 0
      cell_count <- cell_count + 1
      #### call mask strand channel
      strands <- make_strand_mask(offset_px, stage, img_orig, disc_size)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      ### call mask foci foci channel
      bg <- mean(img_orig_foci)
      #### normalise the foci image
      foci_label <- make_foci_mask(offset_factor,bg,crowded_foci,img_orig_foci,brush_size,brush_sigma)
      foci_label <- channel(foci_label, "grey")
      foci_candidates <- computeFeatures.shape(foci_label)
      foci_candidates <- data.frame(foci_candidates)
      foci_areas <- foci_candidates$s.area
      coincident_foci <- get_overlap_mask(strands, foci_label, watershed_stop, img_orig_foci, watershed_radius, watershed_tol)
      foci_per_cell <- get_foci_per_cell(img_file,offset_px,stage,strands,watershed_stop,foci_label,annotation, min_foci,cell_count,img_orig, img_orig_foci, artificial_amp_factor,coincident_foci)
      #### end foci counting function
      if(annotation=="on"){
        cat("\n which counts this many foci:",foci_per_cell, sep = " ")
      }
      image_mat <- as.matrix(img_orig_foci)
      image_mat <- image_mat[image_mat > 1e-06]
      mean_ratio <- median(image_mat)/mean(image_mat)
      skew <- (median(image_mat)-mean(image_mat))/sd(image_mat)
      ### look at properties of the overlap foci.

      ##
      overlap_candidates <- computeFeatures.shape(coincident_foci)
      overlap_candidates <- data.frame(overlap_candidates)
      overlap <- overlap_candidates$s.area
      ### number of foci NOT on an SC
      alone_foci <- nrow(foci_candidates) - foci_per_cell
      if(annotation == "on"){
        cat("\n number of alone foci",alone_foci, sep = " ")
        if(alone_foci < 0){
          cat("\n this one had a negative lone foci amount. Suspect overcounting of foci in this one:")
          plot(rgbImage(strands,foci_label,0*foci_label))
          alone_foci <- 0
        }
      }
      percent_px <- sum(overlap)/sum(foci_areas)
      ### if percent_px too small:
      # count_low_contrast_image()
      tryCatch({
        if(percent_px < 0.3){
          cat("\n file name", file_foci, "coincides", percent_px*100,
              "percent. So decided to reduce the local background thresholding
              requirements.", sep = " ")
          foci_candidates <- make_foci_mask(0.5*offset_factor,0.5*bg,crowded_foci,img_orig_foci,brush_size,brush_sigma)
          ### if it's zero percent... might be genuinely zero. Ignore?
          ### if it's non zero, largish: probably noisy. No distinction between
          ### high background/ high signal vs low bg low signal..
        }
      },
      error = function(e) {
        #what should be done in case of exception?
        cat("\n percent on was a NaN. Suspect no/low signal")
      }
      )
      df_cells <- append_data_frame(WT_str,KO_str,WT_out,KO_out,img_file,foci_areas,df_cells,cell_count,stage,foci_per_cell,image_mat,percent_px,alone_foci)
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
#'
#'
annotate_foci_counting <- function(img_file,cell_count,img_orig,img_orig_foci,artificial_amp_factor,strands,coincident_foci, foci_label){
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
#' @param stage meiotic stage
#' @param foci_per_cell foci count for cell
#' @param image_mat matrix with all pixel values above zero
#' @param percent_px percentage of foci mask that coincides with strand channel
#' small number indicates potentially problematic image.
#' @param alone_foci estimated number of foci that are NOT on a strand.
#'
#'
append_data_frame <- function(WT_str,KO_str,WT_out,KO_out,img_file,foci_areas,df_cells,cell_count,stage,foci_per_cell,image_mat,percent_px,alone_foci){
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
    df_cells <- rbind(df_cells,t(c(img_file,cell_count,genotype,stage,foci_per_cell, sd(foci_areas),mean(foci_areas),median(foci_areas),mean(image_mat),median(image_mat),percent_px,sd(image_mat),alone_foci)))
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
#' @param mask current foci mask
#'
#'
remove_XY <- function(){
  ### if the standard deviation was too big
  ### loop over and delete the giant blobs.
  ### after this function is run, coincident foci needs to be computed again.

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
#'
#'
#'
make_foci_mask <- function(offset_factor,bg,crowded_foci,img_orig_foci,brush_size,brush_sigma){
  foci_mask_crop <- img_orig_foci
  offset <- offset_factor*bg
  if(crowded_foci == TRUE){
    foci_th <- foci_mask_crop > bg + offset
  }
  else{
    ### smooth it
    img_tmp_contrast <- foci_mask_crop
    w <- makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
    img_flo <- filter2(img_tmp_contrast, w)
    ## smooth foci channel
    foci_th <- img_flo > bg + offset
  }
  foci_label <- bwlabel(foci_th)
  return(foci_label)
  # here check whether there is a huge blob. Remove it
  # remove_XY()
  #

}

#' make_strand_mask
#'
#' creates strand mask for strand channel crop
#'
#' @param offset_px, Pixel value offset used in thresholding of dna channel
#' @param stage meitoic stage, currently pachytene or not.
#' @param img_orig original strand crop
#' @param disc_size size of disc for local background calculation in dna channel
#'
#'
make_strand_mask <- function(offset_px, stage, img_orig, disc_size){
  new_img<-img_orig
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
#' @param watershed_stop Turn off default watershed method with "off"
#' @param foci_label black and white mask of foci channel
#' @param img_orig_foci cropped foci channel
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
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
#' @param offset_px, Pixel value offset used in thresholding of dna channel
#' @param stage meitoic stage, currently pachytene or not.
#' @param strands black white mask of strand channel
#' @param watershed_stop Turn off default watershed method with "off"
#' @param foci_label black and white mask of foci channel
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param min_foci minimum pixel area for a foci. Depends on your dpi etc. Defaults to 4
#' @param cell_count unique cell counter
#' @param img_orig original strand crop
#' @param img_orig_foci cropped foci channel
#' @param artificial_amp_factor amplification factor
#' @param coincident_foci mask of coincident foci
#'
#'
#'
get_foci_per_cell <- function(img_file,offset_px,stage,strands,watershed_stop,foci_label, annotation, min_foci, cell_count, img_orig, img_orig_foci, artificial_amp_factor, coincident_foci){
  coincident_df <- data.frame(computeFeatures.shape(coincident_foci))
  if(annotation == "on"){
    print(coincident_df)
  }
  coincident_df <- coincident_df[coincident_df$s.area > min_foci,]
  ### multiply strands by foci_label
  if(annotation == "on"){
    annotate_foci_counting(img_file,cell_count,img_orig,img_orig_foci,artificial_amp_factor,strands,coincident_foci, foci_label)
  }
  foci_per_cell <- nrow(coincident_df)
  return(foci_per_cell)
}
