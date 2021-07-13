#' measure_distances_general
#'
#' Measure the distance between foci on a synaptonemal complex
#'
#' This function first creates masks of synaptonemal complex (SC) and foci
#' channels, where it is assumed that chromosomes are well spread.
#' The user can specify if crowded_foci = TRUE (i.e. >100 foci per cell or so)
#' so that the foci mask is determined with a watershed method.
#' Using the SC mask, individual SCs are isolated and the number of
#' colocalizing foci on a particular strand is determined
#' (in get_distance_general). If a strand meets eccentricity requirements
#' (from EBImage's computeFeatures) and has more than one foci, a middle pixel
#' location is determined with computeFeatures. From here, the function starts
#' counting the distance along the strand in opposite directions
#' (find_start, first_shot_out)
#' It proceeds (get_next_first_dir, get_next_second_dir) by following the
#' maximum intensity in front of it / a gradient
#' approach. These two branches continue until the options to move forward are
#' too dark, i.e. at the end of the SC where a branch will terminate.
#'
#' After a line has been traced out representing an SC, the closest locations
#' of the colocalizing foci to this line are determined. Starting at one end
#' of the line, it starts counting again, recording the pixel distance along
#' that it passes a foci location (get_distances_along).
#'
#' This information is output as a data frame (filename, cell, foci per strand,
#' distances of each along)
#'
#' @export
#' @param img_path, path containing image data to analyse
#' @param stage, meiosis stage of interest. Currently count_foci determines
#' this with thresholding/ object properties in the dna channel. But will be
#' classified using ML model in future versions.
#' @param offset_px, Pixel value offset used in thresholding of dna channel
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be
#' small to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @param foci_norm, Mean intensity to normalise all foci channels to.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param eccentricity_min, The minimum eccentricity (from computefeatures)
#' of a strand to proceed with measuring
#' @param max_strand_area, Maximum pixel area of a strand
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
#' @param watershed_radius Radius (ext variable) in watershed method used in
#'  foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method.
#' Defaults to 0.05.
#' @param crowded_foci TRUE or FALSE, defaults to FALSE.
#' Set to TRUE if you have foci > 100 per cell or so.
#' @param target_foci_number User can specify the number of foci on a strand.
#' Otherwise defaults to any number.
#' @param max_dist_sq maximum distance (squared) between the centre of a foci
#' vs its projected location on the line. Defaults to 20.
#' @param SC_intens_stop mean intensity of a forward moving branch to
#' terminate measuring
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' df_dist <- measure_distances_general(demo_path,offset_factor = 5,
#' brush_size = 1, brush_sigma = 1, annotation = "on", stage = "pachytene")
#' @author Lucy McNeill
#' @return Data frame with properties of synaptonemal (SC) measurements

# should take in same values as count_foci..
measure_distances_general <- function(img_path,offset_px = 0.2, offset_factor = 3, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off", stage = "pachytene", eccentricity_min = 0.6, max_strand_area = 300, channel2_string = "SYCP3", channel1_string = "MLH3",file_ext = "jpeg",KO_str = "--",WT_str = "++",KO_out = "-/-", WT_out = "+/+",watershed_stop = "off", watershed_radius = 1, watershed_tol = 0.05,crowded_foci = FALSE,target_foci_number = FALSE, max_dist_sq = 20, SC_intens_stop = 0.7)
{
  cell_count <- 0
  image_count <-0
  foci_counts <- 0
  foci_count_strand <- c()
  strand_iter <- 0
  dimensionless_dist <- c()
  antibody1_store <- 0
  antibody2_store <- 0
  img_path_new <- paste0(img_path,"/crops/",stage,"/")
  file_list <- list.files(img_path_new)
  df_cols <- c("file","cell_id","genotype","strand_iter","foci_per_strand","iteration_on_strand","foci_x","foci_y","foci_x_line","foci_y_line","total_length", "distance_squared","distance_along","SC_pass_fail")
  df_lengths <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  for (img_file in file_list){
    filename_path_test <- paste0(img_path,"/crops/",stage,"/", img_file)
    img_file <- filename_path_test
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
      file_dna <- img_file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(2*image, "grey")
      antibody1_store <- 1
    }
    if(grepl(paste0('*',channel1_string,'.',file_ext,'$'), img_file)){
      file_foci <- img_file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      antibody2_store <- 1
    }
    if(antibody1_store +antibody2_store == 2){
      antibody1_store <- 0
      antibody2_store <- 0
      cell_count <- cell_count +1
      new_img<-img_orig
      #### now see which have the right amount of strands
      disc <- makeBrush(21, "disc")
      disc <- disc / sum(disc)
      localBackground <- filter2(new_img, disc)
      offset <- offset_px
      thresh_crop <- (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      foci_mask_crop <- img_orig_foci
      bg <- mean(img_orig_foci)
      orig_mean <- mean(img_orig_foci)
      mean_factor <- foci_norm/orig_mean
      img_orig_foci <- img_orig_foci*mean_factor
      foci_label <- threshold_foci_crop(img_orig_foci,offset_factor, brush_size, brush_sigma,stage)
      ##### print properties of the images
      ### multiply strands by foci_label
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      coincident_foci <- bwlabel(foci_label*strands)
      overlap_no <- table(coincident_foci)
      foci_per_cell <-  length(overlap_no)
      image_mat <- as.matrix(foci_mask_crop)
      image_mat <- image_mat[image_mat > 1e-06]
      mean_ratio <- median(image_mat)/mean(image_mat)
      skew <- (median(image_mat)-mean(image_mat))/sd(image_mat)
      ### look at properties of the foci.
      foci_candidates <- computeFeatures.shape(foci_label)
      foci_candidates <- data.frame(foci_candidates)
      foci_areas <- foci_candidates$s.area
      if(annotation == "on"){
        cat("\n looking at cell number:", cell_count, sep = " ")
      }
      dimensionless_dist <- get_distance_general(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,img_file,annotation,eccentricity_min, max_strand_area,cell_count,KO_str ,WT_str,KO_out, WT_out,target_foci_number, max_dist_sq,SC_intens_stop)
      df_lengths <- rbind(df_lengths, dimensionless_dist)
    }
  }
  colnames(df_lengths) <- df_cols
  return(df_lengths)

}

#' get_distance_general
#' Creates mask for SC channel
#'
#' @param strands, A black white mask with SCs as objects
#' @param num_strands, Number of individual strands on SC mask
#' @param new_img, Original strand/dna/SYCP3 channel image with noise removed.
#' @param foci_label, A black white mask with foci as objects
#' @param foci_count_strand, Number of foci counted located on the one SC
#' @param strand_iter, Strand number in iteration over all in cell
#' @param img_file, original filename that cell candidate came from.
#' Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param eccentricity_min, The minimum eccentricity (from computefeatures)
#' of a strand to proceed with measuring
#' @param max_strand_area, Maximum pixel area of a strand
#' @param cell_count Unique cell counter
#' @param KO_str string in filename corresponding to knockout genotype.
#' Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype.
#' Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout.
#' Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout.
#' Defaults to +/+.
#' @param target_foci_number User can specify the number of foci on a strand.
#' Otherwise defaults to any number.
#' @param max_dist_sq maximum distance (squared) between the centre of a foci
#' vs its projected location on the line. Defaults to 20.
#' @param SC_intens_stop mean intensity of a forward moving branch to
#' terminate measuring
#' @return Data frame with properties of synaptonemal (SC) measurements
#'
get_distance_general <- function(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,img_file,annotation, eccentricity_min, max_strand_area,cell_count,KO_str ,WT_str,KO_out, WT_out, target_foci_number, max_dist_sq,SC_intens_stop){
  tryCatch({
    no_strands <- nrow(num_strands)
    strand_count<- 0
    df_cols <- c("file","cell_id","genotype","strand_iter","foci_per_strand","iteration_on_strand","foci_x","foci_y","foci_x_line","foci_y_line","total_length", "distance_squared","distance_along","SC_pass_fail")
    dimensionless_dist <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
    while(strand_count<no_strands){
      strand_iter <- strand_iter + 1
      strand_count <- strand_count + 1
      if (as.numeric(num_strands$s.area[strand_count])<max_strand_area & as.numeric(num_strands$s.area[strand_count])>10){
        tmp_img <- strands
        counter_single <- 0
        # looping over all other objects to crop
        while(counter_single < no_strands){
          counter_single <- counter_single + 1
          # iteratively remove all other objects
          if(counter_single != strand_count){
            tmp_img <- as.matrix(tmp_img)*as.numeric(rmObjects(strands, counter_single, reenumerate = TRUE))
          }
        }
        noise_gone <- bwlabel(tmp_img)*as.matrix(new_img)
        per_strand <- bwlabel(tmp_img)*as.matrix(foci_label)
        per_strand_obj <- computeFeatures.shape(bwlabel(per_strand))
        foci_count_strand <- append(foci_count_strand,nrow(per_strand_obj))
        single_info <- computeFeatures.shape(bwlabel(tmp_img))
        basic_info <- computeFeatures.basic(bwlabel(tmp_img),as.matrix(noise_gone))
        all_info <- computeFeatures(bwlabel(tmp_img),as.matrix(noise_gone))
        moment_info <- computeFeatures.moment(bwlabel(tmp_img),as.matrix(noise_gone))
        moment_info <- as.data.frame(moment_info)
        cx <- moment_info$m.cx
        cy <- moment_info$m.cy
        if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
          genotype <- WT_out
        }
        if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
          genotype <- KO_out
        }
        target_foci_number <- as.integer(target_foci_number)
        if (is.integer(nrow(per_strand_obj))){
          if (target_foci_number > 0){
            if(moment_info$m.eccentricity > eccentricity_min && nrow(per_strand_obj)==target_foci_number){
              ### call the general distance function.
              ## draw box around the middle, find max, locally, use noise_gone as original
              walkers <- 0*noise_gone
              noise_gone <- 2*noise_gone
              window <- 2
              window2 <- window*3
              bright_loc <- find_start(window,noise_gone,cx,cy)
              mean_x <- as.numeric(bright_loc[1,1]) +cx -window2-1
              mean_y <- as.numeric(bright_loc[1,2]) +cy-window2-1
              ix <- (round(mean_x)-window):(round(mean_x)+window)
              iy <- (round(mean_y)-window):(round(mean_y)+window)
              chosen_dir <- get_first_dir(noise_gone,ix,iy,window)
              walkers[round(mean(ix)),round(mean(iy))] <- 1
              ix1 <- ix
              ix2 <- ix
              iy1 <- iy
              iy2 <- iy
              distance_strand <- 0
              distance_strand_2 <- 0
              first_step <- first_shot_out(chosen_dir, ix1,ix2,iy1,iy2,distance_strand, distance_strand_2)
              next_cord <- 2*window+1
              ix1 <- first_step[1:next_cord]
              ix2 <- first_step[(next_cord+1):(2*next_cord)]
              iy1 <- first_step[(2*next_cord+1):(3*next_cord)]
              iy2 <- first_step[(3*next_cord+1):(4*next_cord)]
              distance_strand <- first_step[(4*next_cord+1)]
              distance_strand_2 <- first_step[(4*next_cord+2)]
              dir_1 <- first_step[(4*next_cord+3)]
              dir_2 <- first_step[(4*next_cord+4)]
              new_square_1 <-  noise_gone[ix1,iy1]
              new_square_2 <-  noise_gone[ix2,iy2]
              walkers[round(mean(ix1)),round(mean(iy1))] <- 1
              walkers[round(mean(ix2)),round(mean(iy2))] <- 1
              ## take step in the opposite direction. record new coordinates.
              ## now loop in both directions
              ## set directions to "not done yet" = 0. "done" = 1
              first_dir <- 0
              second_dir <- 0
              start_dir <- chosen_dir
              while(first_dir == 0){
                start_x <- round(mean(ix1))
                start_y <- round(mean(iy1))
                ### call get_first
                dir_1_out <- get_next_first_dir(new_square_1,ix1,iy1,dir_1,window,chosen_dir,distance_strand,first_dir,SC_intens_stop)
                ix1 <- dir_1_out[1:next_cord]
                iy1 <- dir_1_out[(next_cord+1):(2*next_cord)]
                distance_strand <- dir_1_out[(2*next_cord+2)]
                dir_1 <- dir_1_out[(2*next_cord+1)]
                first_dir <- dir_1_out[(2*next_cord+3)]
                start_dir <- dir_1_out[(2*next_cord+4)]
                walkers[round(mean(ix1)),round(mean(iy1))] <- 1
                ## make the new cropped image.
                new_square_1 <-  noise_gone[ix1,iy1]
                if(distance_strand >100){
                  first_dir <- 1
                }
              }

              while(second_dir == 0){
                start_x2 <- round(mean(ix2))
                start_y2 <- round(mean(iy2))
                dir_2_out <- get_next_second_dir(new_square_2,ix2,iy2,dir_2,window,chosen_dir,distance_strand_2,second_dir,SC_intens_stop)
                ix2 <- dir_2_out[1:next_cord]
                iy2 <- dir_2_out[(next_cord+1):(2*next_cord)]
                distance_strand_2 <- dir_2_out[(2*next_cord+2)]
                dir_2 <- dir_2_out[(2*next_cord+1)]
                second_dir <- dir_2_out[(2*next_cord+3)]
                start_dir2 <- dir_2_out[(2*next_cord+4)]
                walkers[round(mean(ix2)),round(mean(iy2))] <- 1
                ## make the new cropped image.
                new_square_2 <-  noise_gone[ix2,iy2]
                if(distance_strand_2 >100){
                  second_dir <- 1
                }
              }
              ### call get_distance_along
              dimensionless_dist <- rbind(dimensionless_dist,get_distances_along(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_count,img_file,annotation,cell_count,strand_iter, per_strand_obj,KO_str ,WT_str,KO_out, WT_out,max_dist_sq))
            }

          }
          else if(moment_info$m.eccentricity > eccentricity_min && nrow(per_strand_obj)>1){
            ## draw box around the middle, find max, locally, use noise_gone as original
            walkers <- 0*noise_gone
            noise_gone <- 2*noise_gone
            window <- 2
            window2 <- window*3
            bright_loc <- find_start(window,noise_gone,cx,cy)
            mean_x <- as.numeric(bright_loc[1,1]) +cx -window2-1
            mean_y <- as.numeric(bright_loc[1,2]) +cy-window2-1
            ix <- (round(mean_x)-window):(round(mean_x)+window)
            iy <- (round(mean_y)-window):(round(mean_y)+window)
            chosen_dir <- get_first_dir(noise_gone,ix,iy,window)
            walkers[round(mean(ix)),round(mean(iy))] <- 1
            ix1 <- ix
            ix2 <- ix
            iy1 <- iy
            iy2 <- iy
            distance_strand <- 0
            distance_strand_2 <- 0
            first_step <- first_shot_out(chosen_dir, ix1,ix2,iy1,iy2,distance_strand, distance_strand_2)
            next_cord <- 2*window+1
            ix1 <- first_step[1:next_cord]
            ix2 <- first_step[(next_cord+1):(2*next_cord)]
            iy1 <- first_step[(2*next_cord+1):(3*next_cord)]
            iy2 <- first_step[(3*next_cord+1):(4*next_cord)]
            distance_strand <- first_step[(4*next_cord+1)]
            distance_strand_2 <- first_step[(4*next_cord+2)]
            dir_1 <- first_step[(4*next_cord+3)]
            dir_2 <- first_step[(4*next_cord+4)]
            new_square_1 <-  noise_gone[ix1,iy1]
            new_square_2 <-  noise_gone[ix2,iy2]
            walkers[round(mean(ix1)),round(mean(iy1))] <- 1
            walkers[round(mean(ix2)),round(mean(iy2))] <- 1
            ## take step in the opposite direction. record new coordinates.
            ## now loop in both directions
            ## set directions to "not done yet" = 0. "done" = 1
            first_dir <- 0
            second_dir <- 0
            start_dir <- chosen_dir
            while(first_dir == 0){
              start_x <- round(mean(ix1))
              start_y <- round(mean(iy1))
              ### call get_first
              dir_1_out <- get_next_first_dir(new_square_1,ix1,iy1,dir_1,window,chosen_dir,distance_strand,first_dir,SC_intens_stop)
              ix1 <- dir_1_out[1:next_cord]
              iy1 <- dir_1_out[(next_cord+1):(2*next_cord)]
              distance_strand <- dir_1_out[(2*next_cord+2)]
              dir_1 <- dir_1_out[(2*next_cord+1)]
              first_dir <- dir_1_out[(2*next_cord+3)]
              start_dir <- dir_1_out[(2*next_cord+4)]
              walkers[round(mean(ix1)),round(mean(iy1))] <- 1
              ## make the new cropped image.
              new_square_1 <-  noise_gone[ix1,iy1]
              if(distance_strand >100){
                first_dir <- 1
              }
            }
            while(second_dir == 0){
              start_x2 <- round(mean(ix2))
              start_y2 <- round(mean(iy2))
              dir_2_out <- get_next_second_dir(new_square_2,ix2,iy2,dir_2,window,chosen_dir,distance_strand_2,second_dir,SC_intens_stop)
              ix2 <- dir_2_out[1:next_cord]
              iy2 <- dir_2_out[(next_cord+1):(2*next_cord)]
              distance_strand_2 <- dir_2_out[(2*next_cord+2)]
              dir_2 <- dir_2_out[(2*next_cord+1)]
              second_dir <- dir_2_out[(2*next_cord+3)]
              start_dir2 <- dir_2_out[(2*next_cord+4)]
              walkers[round(mean(ix2)),round(mean(iy2))] <- 1
              ## make the new cropped image.
              new_square_2 <-  noise_gone[ix2,iy2]
              if(distance_strand_2 >100){
                second_dir <- 1
              }
            }
            ### call get_distance_along
            dimensionless_dist <- rbind(dimensionless_dist,get_distances_along(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_count,img_file,annotation,cell_count,strand_iter, per_strand_obj,KO_str ,WT_str,KO_out, WT_out,max_dist_sq))
          }
          ### then the number of foci was 1
          else if(nrow(per_strand_obj)==1){
            dimensionless_dist <- rbind(dimensionless_dist,c(img_file,cell_count,genotype,strand_iter,1,1,"NA","NA","NA","NA","NA","NA", "NA","NA"))
            colnames(dimensionless_dist) <- df_cols
          }
        }
        ### then the number of foci was zero
        else{
          # a strand with zero foci
          if(target_foci_number == FALSE){
            dimensionless_dist <- rbind(dimensionless_dist,c(img_file,cell_count,genotype,strand_iter,0,"NA","NA","NA","NA","NA","NA","NA", "NA","NA"))
            colnames(dimensionless_dist) <- df_cols
          }
        }
      }
    }
    return(dimensionless_dist)
  },

  error = function(e) {
    str(e) # #prints structure of exception
    if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
      genotype <- WT_out
    }

    if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
      genotype <- KO_out
    }
  }
  )
}

#' get_distances_along
#'
#' Calculates the pixel distance
#'
#' @param distance_strand, total distance along first branch
#' @param distance_strand_2, total distance along second branch
#' @param per_strand, Mask with colocalizing foci
#' @param foci_label, black white mask with foci as objects, not necessarily
#' on SC. Needed for computefeatures.
#' @param walkers, black white mask containing the line that traces through
#' the middle of the SC. Computed earlier in get_distance
#' @param noise_gone, Black white mask with the single SC to be measured
#' @param start_x, x pixel location of where first branch terminated
#' @param start_y, y pixel location of where first branch terminated
#' @param start_x2, x pixel location of where second branch terminated
#' @param start_y2, y pixel location of where second branch terminated
#' @param start_dir, direction that the first branch was traveling along
#' when it terminated
#' @param cx, centre of intensity x of the SC from computefeatures
#' (annotation purposes only)
#' @param cy, centre of intensity y of the SC from computefeatures
#' (annotation purposes only)
#' @param mean_x, starting point x that the two branches move away from
#' to trace out the SC (annotation purposes only)
#' @param mean_y, starting point x that the two branches move away from
#' to trace out the SC (annotation purposes only)
#' @param strand_iter, Strand number in iteration over all in cell
#' @param img_file, original filename that cell candidate came from.
#' Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param cell_count Unique cell number
#' @param uid_strand Unique strand number
#' @param  per_strand_object Foci per strand object
#' @param KO_str string in filename corresponding to knockout genotype.
#' Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype.
#' Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout.
#' Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout.
#' Defaults to +/+.
#' @param max_dist_sq maximum distance (squared) between the centre of a
#' foci vs its projected location on the line. Defaults to 20.
#' @return List of fractional distances between foci for all SCs with two.
#' Optional: total distances of SCs. Optional: images of all resulting traces/ foci locations.
#'
get_distances_along <- function(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_iter,img_file,annotation,cell_count, uid_strand,per_strand_object,KO_str ,WT_str,KO_out, WT_out, max_dist_sq){
  strand_info <- computeFeatures.moment(bwlabel(per_strand),as.matrix(foci_label))
  strand_info <- as.data.frame(strand_info)
  #### turn this into a table/ data frame
  no_foci <- nrow(strand_info)
  cat("\n the number of foci in this strand is", no_foci, sep = " ")
  df_col <- c("file","cell_id","genotype","strand_iter","foci_per_strand","iteration_on_strand","foci_x","foci_y","foci_x_line","foci_y_line","total_length","distance_squared","distance_along","SC_pass_fail")
  foci_df <- data.frame(matrix(ncol = length(df_col), nrow = 0))
  my_walkers_matrix <- t(as.matrix(walkers))
  strand_test <- "pass"
  iter <- 0
  while(iter<= no_foci-1){
    iter <- iter + 1
    cat("\n looking at position data of foci number",iter, sep = " ")
    foci_x <- strand_info$m.cx[iter]
    foci_y <- strand_info$m.cy[iter]
    #### finding closest point
    my_distance_matrix_fi <- 0*my_walkers_matrix+100
    for(row in 1:nrow(my_walkers_matrix)) {
      for(col in 1:ncol(my_walkers_matrix)) {
        if(my_walkers_matrix[row, col]==1){
          my_distance_matrix_fi[row, col] <- (row-strand_info$m.cy[iter])^2+(col-strand_info$m.cx[iter])^2
        }
      }
    }
    bright_loc_fi <- which(my_distance_matrix_fi == min(my_distance_matrix_fi),arr.ind = TRUE)
    mean_x <- as.numeric(bright_loc_fi[1,1])
    mean_y <- as.numeric(bright_loc_fi[1,2])
    distance_fi <- (strand_info$m.cy[iter]-mean_x)^2 +(strand_info$m.cx[iter]-mean_y)^2
    if(distance_fi > max_dist_sq){
      strand_test <- "fail"
    }
    if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
      genotype <- WT_out
    }
    if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
      genotype <- KO_out
    }
    foci_df <- rbind(foci_df, c(img_file,cell_count,genotype,strand_iter,no_foci,iter,foci_x,foci_y,mean_y,mean_x,(distance_strand+distance_strand_2),distance_fi, "NA","NA"))
    ##### now count dimensionless distance
    x_curr <- start_x
    y_curr <- start_y
    dir_curr <- start_dir
    if (dir_curr <5){
      dir_curr <- dir_curr +4
    }
    else{
      dir_curr <- dir_curr -4
    }
    looping <- 1
    test_walker <-  0*my_walkers_matrix
    length_walker <- 1
    ## measuring foci
    dist_along <- 0
    while(looping == 1){
      length_walker <- length_walker + 1
      if (length_walker > 150){
        print("probably an infinite loop")
        strand_test <- "fail"
        looping <- 0
      }
      #### make the new df row with distance along here.
      colnames(foci_df) <- df_col
      if(x_curr == as.numeric(foci_df$foci_x_line[iter]) && y_curr == as.numeric(foci_df$foci_y_line[iter])){
        ### record the distance along
        foci_dist_along <- dist_along
        cat("\n I found a focus, which is this many pixels along:", dist_along, sep = " ")
        cat("\n and its distance (squared) to the actual foci centre was", distance_fi, sep = " ")
        cat("\n and this was a pass/fail", strand_test,  sep = " ")
        if(annotation == "on"){
          ch1 <- bwlabel(noise_gone)
          ch1 <-channel(noise_gone,"grey")
          ch2 <- bwlabel(foci_label)
          ch2 <- channel(foci_label,"grey")
          ch3 <- channel(bwlabel(walkers),"grey")
          bluered <- rgbImage(ch1, ch2, 0*ch1)
          plot(bluered)
          print("with foci locations included (actual)")
          plot(bluered)
          text(x = foci_x, y = foci_y, label = "+", col = "blue", cex = 2)
          print("with foci locations included, on the line")
          plot(bluered)
          text(x = mean_y, y = mean_x, label = "+", col = "magenta", cex = 2)
          blueredgreen <- rgbImage(ch1, ch2, ch3)
          plot(blueredgreen)
          text(x = mean_y, y = mean_x, label = "+", col = "magenta", cex = 2)
        }
        looping <- 0
        foci_df$distance_along[iter] <- dist_along
        foci_df$genotype[iter] <- genotype
        foci_df$SC_pass_fail[iter] <- strand_test
      }
      ### do the opposite of the travelling direction
      if(dir_curr == 1){
        # 8 position
        if(walkers[x_curr-1,y_curr-1]==1){
          y_curr <- y_curr -1
          x_curr <- x_curr -1
          dir_curr <- 8
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 1 position
        else if(walkers[x_curr,y_curr-1]==1){
          y_curr <- y_curr -1
          dir_curr <- 1
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 2 position
        else if(walkers[x_curr+1,y_curr-1]==1){
          y_curr <- y_curr -1
          x_curr <- x_curr +1
          dir_curr <- 2
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
      }
      else if(dir_curr == 2){
        # 1 position
        if(walkers[x_curr,y_curr-1]==1){
          y_curr <- y_curr -1
          dir_curr <- 1
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 2 position
        else if(walkers[x_curr+1,y_curr-1]==1){
          y_curr <- y_curr -1
          x_curr <- x_curr +1
          dir_curr <- 2
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 3 position
        else if(walkers[x_curr+1,y_curr]==1){
          x_curr <- x_curr +1
          dir_curr <- 3
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
      }
      else if(dir_curr == 3){
        # 2 position
        if(walkers[x_curr+1,y_curr-1]==1){
          y_curr <- y_curr -1
          x_curr <- x_curr +1
          dir_curr <- 2
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 3 position
        else if(walkers[x_curr+1,y_curr]==1){
          x_curr <- x_curr +1
          dir_curr <- 3
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 4 position
        else if(walkers[x_curr+1,y_curr+1]==1){
          y_curr <- y_curr +1
          x_curr <- x_curr +1
          dir_curr <- 4
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
      }
      else if(dir_curr == 4){
        # 3 position
        if(walkers[x_curr+1,y_curr]==1){
          x_curr <- x_curr +1
          dir_curr <- 3
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 4 position
        else if(walkers[x_curr+1,y_curr+1]==1){
          y_curr <- y_curr +1
          x_curr <- x_curr +1
          dir_curr <- 4
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 5 position
        else if(walkers[x_curr,y_curr+1]==1){
          y_curr <- y_curr +1
          dir_curr <- 5
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
      }
      else if(dir_curr == 5){
        # 4 position
        if(walkers[x_curr+1,y_curr+1]==1){
          y_curr <- y_curr +1
          x_curr <- x_curr +1
          dir_curr <- 4
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 5 position
        else if(walkers[x_curr,y_curr+1]==1){
          y_curr <- y_curr +1
          dir_curr <- 5
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 6 position
        else if(walkers[x_curr-1,y_curr+1]==1){
          x_curr <- x_curr -1
          y_curr <- y_curr +1
          dir_curr <- 6
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
      }
      else if(dir_curr == 6){
        # 5 position
        if(walkers[x_curr,y_curr+1]==1){
          y_curr <- y_curr +1
          dir_curr <- 5
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 6 position
        else if(walkers[x_curr-1,y_curr+1]==1){
          x_curr <- x_curr -1
          y_curr <- y_curr +1
          dir_curr <- 6
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 7 position
        else if(walkers[x_curr-1,y_curr]==1){
          x_curr <- x_curr -1
          dir_curr <- 7
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
      }
      else if(dir_curr == 7){
        # 6 position
        if(walkers[x_curr-1,y_curr+1]==1){
          x_curr <- x_curr -1
          y_curr <- y_curr +1
          dir_curr <- 6
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 7 position
        else if(walkers[x_curr-1,y_curr]==1){
          x_curr <- x_curr -1
          dir_curr <- 7
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 8 position
        else if(walkers[x_curr-1,y_curr-1]==1){
          x_curr <- x_curr -1
          y_curr <- y_curr -1
          dir_curr <- 8
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
      }
      else if(dir_curr == 8){
        # 7 position
        if(walkers[x_curr-1,y_curr]==1){
          x_curr <- x_curr -1
          dir_curr <- 7
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
        # 8 position
        else if(walkers[x_curr-1,y_curr-1]==1){
          x_curr <- x_curr -1
          y_curr <- y_curr -1
          dir_curr <- 8
          length_walker <- length_walker + 1
          dist_along <- dist_along + sqrt(2)
        }
        # 1 position
        else if(walkers[x_curr,y_curr-1]==1){
          y_curr <- y_curr -1
          dir_curr <- 1
          length_walker <- length_walker + 1
          dist_along <- dist_along + 1
        }
      }
    }
  }
  colnames(foci_df) <- df_col
  if(length_walker<149){
    uid <- strand_iter
    ### add the new row here...
    if(strand_test == "fail"){
      foci_loop <- 1
      while(foci_loop <= no_foci){
        foci_df$SC_pass_fail[foci_loop] <- "fail"
        foci_loop <- foci_loop +1
      }
    }
    return(foci_df)
  }
}

#' threshold_SC_crop
#'
#' Creates mask for SC channel
#'
#' @param image SC channel image
#' @param offset offset for an SC signal
#' @return A black white mask with SCs as objects
#'
threshold_SC_crop <- function(image, offset){
  disc <- makeBrush(21, "disc")
  disc <- disc / sum(disc)
  localBackground <- filter2(image, disc)
  thresh_crop <- (image - localBackground > offset)
  strands <- bwlabel(thresh_crop)
  return(strands)
}


#' threshold_foci_crop
#'
#' Creates mask for foci channel
#'
#' @param image foci channel image
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small
#'  to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be
#' small to avoid erasing foci.
#' @param crowded_foci TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have foci > 100 or so.
#' @return A black white mask with foci as objects
#'
threshold_foci_crop <- function(image, offset_factor, brush_size, brush_sigma, crowded_foci){
  bg <- mean(image)
  offset <- offset_factor*bg

  ### new stuff July
  if(crowded_foci == TRUE){
    foci_th <- image > bg + offset
    #foci_th <- watershed(bwlabel(foci_th)*as.matrix(img_orig_foci),tolerance=0.05, ext=1)
  }
  else{
    ### smooth it
    img_tmp_contrast <- image
    w <- makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
    #w = makeBrush(size = 1, shape = 'gaussian', sigma = 3)
    img_flo <- filter2(img_tmp_contrast, w)
    ## smooth foci channel
    foci_th <- img_flo > bg + offset
  }
  foci_label <- bwlabel(foci_th)
  foci_label <- channel(foci_label, "grey")
  return(foci_label)
}


#' find_start
#'
#' Finds an appropriate place to start counting along an SC by looking in a window around the "centre of intensity"
#'
#' @param window window size that we want to look in
#' @param noise_gone Image containing only the SC of interest without background
#' @param cx,cy The location of the "centre of intensity" using computeFeatures
#' (not necessarily on SC)
#' @return The location of the starting point (on the SC)
#'
find_start <- function(window,noise_gone,cx,cy){
  window2 <- window*3
  ix <- (round(cx)-window):(round(cx)+window)
  iy <- (round(cy)-window):(round(cy)+window)
  ix_c <- (round(cx)-window2):(round(cx)+window2)
  iy_c <- (round(cy)-window2):(round(cy)+window2)
  sub_square <- noise_gone[ix_c,iy_c]
  ### not sure about this
  bright_loc <- which(as.matrix(sub_square) == max(as.matrix(sub_square)),arr.ind = TRUE)
  return(bright_loc)
}


#' get_first_dir
#'
#' Finds an appropriate place to start counting along an SC by looking in a window around the "centre of intensity"
#'
#' @param window window size where we compute local gradients
#' @param noise_gone Image containing only the SC of interest without background
#' @param ix,iy Pixel locations of the starting point (on the SC)
#' @return string: The direction of brightest intensity from starting point
#'
get_first_dir <- function(noise_gone,ix,iy,window){
  sub_square2 <- noise_gone[ix,iy]
  ## shoot out vertical and horizontal lines
  # vertical
  line_1 <- sub_square2[window,]
  # horizontal
  line_2 <- sub_square2[,window]
  ## shoot out diagonals, accounting for factor sqrt(2)
  # clockwise
  line_3 <- diag(sub_square2)
  # anti-clockwise
  line_4 <- diag(rotate(sub_square2, 90))
  max_mean <- max(mean(line_1),mean(line_2),mean(line_3),mean(line_4))
  if(abs(mean(line_2)-max_mean) <1e-4){
    chosen_dir <- "horizontal"
  }

  if(abs(mean(line_1)-max_mean) <1e-4){
    chosen_dir <- "vertical"
  }

  if(abs(mean(line_3)-max_mean) <1e-4){
    chosen_dir <- "diag"
  }
  if(abs(mean(line_4)-max_mean) <1e-4){
    chosen_dir <- "off-diag"
  }
  return(chosen_dir)
}


#' first_shot_out
#'
#' Moves one pixel away from the starting point
#'
#' @param chosen_dir string: brightest direction of a line passing through
#' starting point
#' @param ix1, starting point x
#' @param ix2, starting point x
#' @param iy1, starting point y
#' @param iy2, starting point y
#' @param distance_strand, distance along first branch (zero).
#' @param distance_strand_2, distance along second branch (zero).
#' @return New sub square for first and second branch

first_shot_out <- function(chosen_dir, ix1,ix2,iy1,iy2,distance_strand, distance_strand_2){
  ### take step in one side of the direction chosen from mean.
  if(chosen_dir == "horizontal"){
    ### right branch: move +1 in x
    ## might need to use something different to ix and iy? or start calling them something else?
    dir_1 <- 3
    dir_2 <- 7
    ix1 <- ix1 +1
    ## left branch: move -1 in x
    ix2 <- ix2 -1
    distance_strand <- 1
    distance_strand_2 <- 1
  }

  if(chosen_dir == "vertical"){
    dir_1 <- 1
    dir_2 <- 5
    iy1 <- iy1 +1
    ## left branch: move -1 in x
    iy2 <- iy2 -1
    distance_strand <- 1
    distance_strand_2 <- 1
  }
  if(chosen_dir == "diag"){
    dir_1 <- 4
    dir_2 <- 8
    ix1 <- ix1 +1
    iy1 <- iy1 +1

    distance_strand <- sqrt(2)
    distance_strand_2 <- sqrt(2)

    ## left branch: move -1 in x
    ix2 <- ix2 -1
    iy2 <- iy2 -1
  }

  if(chosen_dir == "off-diag"){
    dir_1 <- 2
    dir_2 <- 6
    ix1 <- ix1 +1
    iy1 <- iy1 -1

    ## left branch: move -1 in x
    ix2 <- ix2 -1
    iy2 <- iy2+1
    distance_strand <- sqrt(2)
    distance_strand_2 <- sqrt(2)
  }
  return(c(ix1,ix2,iy1,iy2,distance_strand,distance_strand_2,dir_1,dir_2))
}




#' get_next_first_dir
#'
#' Moves one pixel away one first branch
#'

#' @param new_square_1, the subsquare that contains the x y location of
#' first branch walker in the middle.
#' @param ix1, current x position of first branch walker along SC
#' @param iy1, current y position of first branch walker along SC
#' @param dir_1, The direction (choice of 8) of the first branch step
#' @param window, number of pixels ahead that the intensity gradient is
#' computed with
#' @param chosen_dir The brightest direction (choice of 4) of the previous step
#' @param distance_strand, current distance along the first branch
#' @param first_dir, zero while still measuring along first branch, one when
#' first branch terminates and counting stops.
#' @param SC_intens_stop mean intensity of a forward moving branch to
#' terminate measuring
#' @return New sub square for first branch after taking one step
#'
get_next_first_dir <- function(new_square_1,ix1,iy1,dir_1,window,chosen_dir,distance_strand,first_dir,SC_intens_stop){
  start_x <- round(mean(ix1))
  start_y <- round(mean(iy1))
  start_dir <- dir_1
  #start_dir2 <- dir_2

  ## make a new crop with the new choice at centre.
  ### looking at new_square_1
  ### you have the direction that you moved in
  ## now dealing with new square (already made at end of last loop...)

  ### now make all the 8 choices
  if (dir_1 ==1){
    ## take the latest crop.
    ## crop it depending on direction
    # only compute mean of one direction..
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(diag(new_square_1[ix_8,iy_8]))
    ## 1 position
    iy_1 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(new_square_1[window+1,iy_1])
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(diag(rotate(new_square_1[ix_2,iy_2], 90)))
    ## determine the maximum mean and corresponding direction
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }

    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "up"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right-up"
    }
  }

  if (dir_1 ==2 ){
    ## 1 position
    iy_1 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(new_square_1[window+1,iy_1])
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(diag(rotate(new_square_1[ix_2,iy_2], 90)))
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(new_square_1[ix_3,window+1])
    ## determine the maximum mean and corresponding direction
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "up"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right-up"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right"
    }
  }
  if (dir_1 == 3 ){
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(diag(rotate(new_square_1[ix_2,iy_2], 90)))
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(new_square_1[ix_3,window+1])
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(diag(new_square_1[ix_4,iy_4]))
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right-up"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }
  }
  if (dir_1 == 4 ){
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(new_square_1[ix_3,window+1])
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(diag(new_square_1[ix_4,iy_4]))
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(new_square_1[window+1,iy_5])
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "down"
    }
  }
  if (dir_1 == 5){
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(diag(new_square_1[ix_4,iy_4]))
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(new_square_1[window+1,iy_5])
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(diag(rotate(new_square_1[ix_6,iy_6], 90)))
    max_mean <- max(mean_1,mean_2,mean_3)

    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }

    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "down"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
  }

  if (dir_1 == 6){
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(new_square_1[window+1,iy_5])
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(diag(rotate(new_square_1[ix_6,iy_6], 90)))
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(new_square_1[ix_7,window+1])
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "down"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left"
    }
  }
  if (dir_1 == 7){
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(diag(rotate(new_square_1[ix_6,iy_6], 90)))
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(new_square_1[ix_7,window+1])
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(diag(new_square_1[ix_8,ix_8]))
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }
  }
  if (dir_1 == 8){
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(new_square_1[ix_7,window+1])
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(diag(new_square_1[ix_8,iy_8]))
    ## 1 position
    iy_1 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(new_square_1[window+1,iy_1])
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "up"
    }
  }
  ## need some condition to stop measuring? If mean is < 0.1? something like this..
  if(max_mean < SC_intens_stop){
    first_dir <- 1
  }
  ## now see which mean/direction is maximised. Choose this
  if(chosen_dir == "up"){
    ## move up in y
    iy1 <- iy1-1
    distance_strand <- distance_strand + 1
    dir_1 <- 1
  }
  if(chosen_dir == "right"){
    ## move up in x
    ix1 <- ix1+1
    distance_strand <- distance_strand + 1
    dir_1 <- 3
  }
  if(chosen_dir == "down"){
    ## move down in y
    iy1 <- iy1+1
    distance_strand <- distance_strand + 1
    dir_1 <- 5
  }
  if(chosen_dir == "left"){
    ## move down in x
    ix1 <- ix1-1
    distance_strand <- distance_strand + 1
    dir_1 <- 7
  }
  ## now see which mean/direction is maximised. Choose this
  if(chosen_dir == "right-up"){
    ## move up in y
    ## move up in x
    iy1 <- iy1-1
    ix1 <- ix1+1
    distance_strand <- distance_strand + sqrt(2)
    dir_1 <- 2
  }
  if(chosen_dir == "right-down"){
    ## move up in x
    ## move down in y
    ix1 <- ix1+1
    iy1 <- iy1+1
    distance_strand <- distance_strand + sqrt(2)
    dir_1 <- 4
  }
  if(chosen_dir == "left-down"){
    ## move down in x
    ## move down in y
    iy1 <- iy1+1
    ix1 <- ix1 -1
    distance_strand <- distance_strand + sqrt(2)
    dir_1 <- 6
  }
  if(chosen_dir == "left-up"){
    ## move down in x
    ## move up in y
    ix1 <- ix1-1
    iy1 <- iy1-1
    distance_strand <- distance_strand + sqrt(2)
    dir_1 <- 8
  }
  return(c(ix1,iy1,dir_1,distance_strand,first_dir,start_dir))

}



#' get_next_second_dir
#'
#' Moves one pixel away one second branch. Terminates if at the end of the SC.
#' @param new_square_2, the subsquare that contains the x y location of second
#' branch walker in the middle.
#' @param ix2, current x position of second branch walker along SC
#' @param iy2, current y position of second branch walker along SC
#' @param dir_2, The direction (choice of 8) of the second (opposite to first)
#' branch step
#' @param window, number of pixels ahead that the intensity gradient is computed
#' with
#' @param chosen_dir The brightest direction (choice of 4) of the previous step
#' @param distance_strand_2, current distance along the second branch
#' @param second_dir, zero while still measuring along first branch, one when
#' first branch terminates and counting stops.
#' @param SC_intens_stop mean intensity of a forward moving branch to terminate
#' measuring

#' @return New sub square for second branch after taking one step
#'

get_next_second_dir <- function(new_square_2,ix2,iy2,dir_2,window,chosen_dir,distance_strand_2,second_dir,SC_intens_stop){
  start_dir2 <- dir_2
  start_x2 <- round(mean(ix2))
  start_y2 <- round(mean(iy2))
  if (dir_2 ==1){
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(diag(new_square_2[ix_8,iy_8]))
    ## 1 position
    iy_1 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(new_square_2[window+1,iy_1])
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(diag(rotate(new_square_2[ix_2,iy_2], 90)))
    ## determine the maximum mean and corresponding direction
    max_mean <- max(mean_1,mean_2,mean_3)

    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "up"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right-up"

    }
  }

  if (dir_2 ==2 ){
    mean_1 <- mean(new_square_2[window,1:window])
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(diag(rotate(new_square_2[ix_2,iy_2], 90)))
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(new_square_2[ix_3,window+1])
    ## determine the maximum mean and corresponding direction
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "up"
    }

    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right-up"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right"
    }

  }

  if (dir_2 == 3 ){
    ## 2 position
    ix_2 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_2 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(diag(rotate(new_square_2[ix_2,iy_2], 90)))
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(new_square_2[ix_3,window+1])
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(diag(new_square_2[ix_4,iy_4]))
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right-up"
    }

    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }
  }

  if (dir_2 == 4 ){
    ## 3 position
    ix_3 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(new_square_2[ix_3,window+1])
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(diag(new_square_2[ix_4,iy_4]))
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(new_square_2[window+1,iy_5])

    max_mean <- max(mean_1,mean_2,mean_3)

    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right"
    }

    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "down"
    }
  }

  if (dir_2 == 5){
    ## 4 position, new one is now in corner, mean of diagonal
    ix_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    iy_4 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(diag(new_square_2[ix_4,iy_4]))
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(new_square_2[window+1,iy_5])
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_3 <- mean(diag(rotate(new_square_2[ix_6,iy_6], 90)))
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "right-down"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "down"
    }

    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
  }
  if (dir_2 == 6){
    ## 5 position
    iy_5 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(new_square_2[window+1,iy_5])
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_2 <- mean(diag(rotate(new_square_2[ix_6,iy_6], 90)))
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(new_square_2[ix_7,window+1])
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "down"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left"
    }
  }
  if (dir_2 == 7){
    ## 6 position
    ix_6 <- seq(from = 1, to = window+1, by = 1)
    iy_6 <- seq(from = window+1, to = 2*window+1, by = 1)
    mean_1 <- mean(diag(rotate(new_square_2[ix_6,iy_6], 90)))
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(new_square_2[ix_7,window+1])
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(diag(new_square_2[ix_8,iy_8]))
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left-down"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }
  }
  if (dir_2 == 8){
    # 7 position
    ix_7 <- seq(from = 1, to = window+1, by = 1)
    mean_1 <- mean(new_square_2[ix_7,window+1])
    ## 8 position
    ix_8 <- seq(from = 1, to = window+1, by = 1)
    iy_8 <- seq(from = 1, to = window+1, by = 1)
    mean_2 <- mean(diag(new_square_2[ix_8,iy_8]))
    ## 1 position
    iy_1 <- seq(from = 1, to = window+1, by = 1)
    mean_3 <- mean(new_square_2[window+1,iy_1])
    max_mean <- max(mean_1,mean_2,mean_3)
    if(abs(mean_1-max_mean) <1e-4){
      chosen_dir <- "left"
    }
    if(abs(mean_2-max_mean) <1e-4){
      chosen_dir <- "left-up"
    }
    if(abs(mean_3-max_mean) <1e-4){
      chosen_dir <- "up"
    }
  }
  ## need condition to stop measuring. SC_intens_stop
  if(max_mean < SC_intens_stop){
    second_dir <- 1
  }
  ## now see which mean/direction is maximised. Choose this
  if(chosen_dir == "up"){
    ## move up in y
    iy2 <- iy2-1
    distance_strand_2 <- distance_strand_2 + 1
    dir_2 <- 1
  }
  if(chosen_dir == "right"){
    ## move up in x
    ix2 <- ix2+1
    distance_strand_2 <- distance_strand_2 + 1
    dir_2 <- 3
  }
  if(chosen_dir == "down"){
    ## move down in y
    iy2 <- iy2+1
    distance_strand_2 <- distance_strand_2 + 1
    dir_2 <- 5
  }
  if(chosen_dir == "left"){
    ## move down in x
    ix2 <- ix2-1
    distance_strand_2 <- distance_strand_2 + 1
    dir_2 <- 7
  }
  ## now see which mean/direction is maximised. Choose this
  if(chosen_dir == "right-up"){
    ## move up in y
    ## move up in x
    iy2 <- iy2-1
    ix2 <- ix2+1
    distance_strand_2 <- distance_strand_2 + sqrt(2)
    dir_2 <- 2
  }
  if(chosen_dir == "right-down"){
    ## move up in x
    ## move down in y
    ix2 <- ix2+1
    iy2 <- iy2+1
    distance_strand_2 <- distance_strand_2 + sqrt(2)
    dir_2 <- 4
  }
  if(chosen_dir == "left-down"){
    ## move down in x
    ## move down in y
    iy2 <- iy2+1
    ix2 <- ix2 -1
    distance_strand_2 <- distance_strand_2 + sqrt(2)
    dir_2 <- 6
  }
  if(chosen_dir == "left-up"){
    ## move down in x
    ## move up in y
    ix2 <- ix2-1
    iy2 <- iy2-1
    distance_strand_2 <- distance_strand_2 + sqrt(2)
    dir_2 <- 8
  }
  return(c(ix2,iy2,dir_2,distance_strand_2,second_dir,start_dir2))
}
