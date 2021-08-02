#' get_pachytene
#'
#' Identifies crops in pachytene
#'
#' This function takes the crops make by auto_crop fast, and determines the
#' number of synaptonemal complex candidates by considering the local background
#' and using EBImage functions. In general, very bright objects which contrast
#' highly with the background will be classified as the same object. Dim objects
#' will likely be classified as many different objects. If the number of objects
#' is too high compared to the species number (species_num) then the cell is
#' determined to not be in pachytene. Note that this function has been optimized
#' for mouse cells which can be very well spread / separated.
#'
#' @export get_pachytene
#' @param img_path, path containing crops analyse
#' @param path_out, user specified output path. Defaults to img_path
#' @param species_num, number of chromosomes in the species
#' @param offset, Pixel value offset used in therholding for the dna (SYCP3) channel
#' @param ecc_thresh, The minimum average eccentricity of all objects in mask determined by computefeatures, for a cell to be pachytene.
#' @param area_thresh, The minimum ratio of pixels included in mask to total, for a cell to be classified as pachytene.
#' @param channel1_string String appended to the files showing the channel illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel illuminating synaptonemal complexes. Defaults to SYCP3
#' @param file_ext file extension of your images e.g. tiff jpeg or png.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param KO_str string in filename corresponding to knockout genotype. Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype. Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout. Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout. Defaults to +/+.
#' @param resize_l length of resized square cell image.
#' @param artificial_amp_factor Amplification of foci channel, for RGB output files. Deaults to 3.
#' @param strand_amp multiplication of strand channel.
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' SYCP3_stats <- get_pachytene(demo_path,ecc_thresh = 0.8, area_thresh = 0.04, annotation = "on")
#' @author Lucy McNeill
#' @return Pairs of foci and SC channel crops for pachytene
#'


get_pachytene <- function(img_path, species_num = 20, offset = 0.2,ecc_thresh = 0.85, area_thresh = 0.06, annotation = "off", channel2_string = "SYCP3", channel1_string = "MLH3",file_ext = "jpeg", KO_str = "--",WT_str = "++",KO_out = "-/-", WT_out = "+/+", path_out = img_path, artificial_amp_factor=3, strand_amp = 2, resize_l = 120)
{
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  pachytene_count <- 0
  max_obj <- species_num + round(0.1*species_num)
  df_cols <- c("filename","cell_no","genotype", "px_mask","px_total", "px_fraction", "mean_ecc","mean_ratio","skew","sd_bright_px","stage_classification")
  df_cells <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_cells) <- df_cols
  img_path_new <- paste0(img_path,"/crops")
  if(img_path == path_out){
    img_path_out <- img_path_new
    dir.create(paste0(img_path_new,"/pachytene"))
    dir.create(paste0(img_path_new,"/pachytene-RGB"))
  }
  else{
    img_path_out <- path_out
    dir.create(paste0(path_out,"/pachytene"))
    dir.create(paste0(path_out,"/pachytene-RGB"))
  }

  file_list <- list.files(img_path_new)
  for (img_file in file_list){
    file_base <- img_file
    filename_path_test <- paste0(img_path,"/crops/", img_file)
    img_file <- filename_path_test
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
      file_dna <- img_file
      file_base_dna <- file_base
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig_highres <- channel(image, "grey")
      img_orig <- channel(strand_amp*image, "grey")
      img_orig <- resize(img_orig, w = resize_l, h = resize_l)
      antibody1_store <- 1
    }
    if(grepl(paste0('*',channel1_string,'.',file_ext,'$'), img_file)){
      file_base_foci <- file_base
      file_foci <- img_file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      antibody2_store <- 1
    }
    if(antibody1_store + antibody2_store ==2){
      antibody1_store <- 0
      antibody2_store <- 0
      new_img<-img_orig
      #### now see which have the right amount of strands
      disc <- makeBrush(21, "disc")
      disc <- disc / sum(disc)
      localBackground <- filter2(new_img, disc)
      thresh_crop <- (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      #### segment the strands
      if (nrow(num_strands)<max_obj && nrow(num_strands)>3){
        cell_count <- cell_count + 1
        ### identified a good image. count foci
        ### data frame stuff
        if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
          genotype <- WT_out
        }
        if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
          genotype <- KO_out
        }
        image_mat <- as.matrix(new_img)
        image_mat <- image_mat[image_mat > 1e-01]
        mean_ratio <- median(image_mat)/mean(image_mat)
        skew <- (median(image_mat)-mean(image_mat))/sd(image_mat)
        ### look at properties of the foci.
        SC_candidates <- computeFeatures.shape(strands)
        SC_candidates <- data.frame(SC_candidates)
        SC_areas <- SC_candidates$s.area
        moment_info <- computeFeatures.moment(bwlabel(strands),as.matrix(new_img))
        moment_info <- as.data.frame(moment_info)
        ## might actually want to find the real centre first..
        ### data frame stuff ends
        dim_img <- nrow(new_img)
        num_px <- dim_img*dim_img
        px_mask <- sum(SC_areas)
        px_total <- num_px
        px_fraction <- px_mask/px_total
        mean_ecc <- mean(moment_info$m.eccentricity)
        stage_classification <- "pachytene"
        sd_bright_px <- sd(image_mat)
        if(mean_ecc > ecc_thresh){
          if(px_fraction > area_thresh){
            stage_classification <- "pachytene"
            df_cells <- rbind(df_cells,t(c(img_file,cell_count,genotype,px_mask, px_total,px_fraction, mean_ecc,mean_ratio,skew,sd_bright_px,stage_classification)))
            pachytene_count <- pachytene_count + 1
            file_dna <- tools::file_path_sans_ext(file_base_dna)
            filename_crop <- paste0(img_path_out,"/pachytene/", file_dna,'.',file_ext)
            writeImage(img_orig_highres, filename_crop)
            if(annotation == "on"){
              print("decided the following is pachytene")
              display(img_orig_highres)
            }
            file_foci <- tools::file_path_sans_ext(file_base_foci)
            filename_crop_foci <- paste0(img_path_out,"/pachytene/", file_foci,'.',file_ext)
            writeImage(img_orig_foci, filename_crop_foci)
            ### add RGB channel
            ch1 <-channel(img_orig_highres,"grey")
            ch2 <- channel(img_orig_foci,"grey")
            RGB_img <- rgbImage(ch1,ch2,0*ch1)
            filename_crop_RGB <- paste0(img_path_out,"/pachytene-RGB/", file_dna,'.',file_ext)
            writeImage(RGB_img, filename_crop_RGB)
          }
        }
      }
    }
  }
cat("number of cells kept",pachytene_count, sep = " ")
colnames(df_cells) <- df_cols
return(df_cells)
}

