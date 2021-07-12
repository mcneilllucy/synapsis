#' measure_distances_general
#'
#' Measure the distance between foci on a synaptonemal complex
#'
#' @export
#' @param img_path, path containing image data to analyse
#' @param stage, meiosis stage of interest. Currently count_foci determines this with thresholding/ object properties in the dna channel. But will be classified using ML model in future versions.
#' @param offset_px, Pixel value offset used in thresholding of dna channel
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be small to avoid erasing foci.
#' @param foci_norm, Mean intensity to normalise all foci channels to.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param eccentricity_min, The minimum eccentricity (from computefeatures) of a strand to proceed with measuring
#' @param max_strand_area, Maximum pixel area of a strand
#' @param channel1_string String appended to the files showing the channel illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel illuminating synaptonemal complexes. Defaults to SYCP3
#' @param file_ext file extension of your images e.g. tiff jpeg or png.
#' @param KO_str string in filename corresponding to knockout genotype. Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype. Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout. Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout. Defaults to +/+.
#' @param watershed_stop Turn off default watershed method with "off"
#' @param watershed_radius Radius (ext variable) in watershed method used in foci channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @return Data frame with properties of synaptonemal (SC) measurements

# should take in same values as count_foci..
measure_distances_general <- function(img_path,offset_px = 0.2, offset_factor = 3, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off", stage = "pachytene", eccentricity_min = 0.6, max_strand_area = 300, channel2_string = "SYCP3", channel1_string = "MLH3",file_ext = "jpeg",KO_str = "--",WT_str = "++",KO_out = "-/-", WT_out = "+/+",watershed_stop = "off", watershed_radius = 1, watershed_tol = 0.05)
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
  #colnames(df_lengths) <- df_cols
  ## for each image that is *-dna.jpeg,
  for (img_file in file_list){
    filename_path_test <- paste0(img_path,"/crops/",stage,"/", img_file)
    img_file <- filename_path_test
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
    #if(grepl("*SYCP3.jpeg", file)){
      file_dna <- img_file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(2*image, "grey")
      antibody1_store <- 1
    }
    if(grepl(paste0('*',channel1_string,'.',file_ext,'$'), img_file)){
    #if(grepl("*MLH3.jpeg", file)){
      file_foci <- img_file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(  antibody1_store +antibody2_store == 2){
      antibody1_store <- 0
      antibody2_store <- 0
      cell_count <- cell_count +1
      new_img<-img_orig
      #display(new_img)
      #### now see which have the right amount of strands
      disc <- makeBrush(21, "disc")
      disc <- disc / sum(disc)
      localBackground <- filter2(new_img, disc)
      offset <- offset_px
      thresh_crop <- (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      #display(strands)
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
        print("looking at cell number:")
        print(cell_count)
      }
      ################ distance starts (make function later)
      dimensionless_dist <- get_distance_general(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,img_file,annotation,eccentricity_min, max_strand_area,cell_count,KO_str ,WT_str,KO_out, WT_out)
      #colnames(dimensionless_dist) <- df_cols
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
#' @param img_file, original filename that cell candidate came from. Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param eccentricity_min, The minimum eccentricity (from computefeatures) of a strand to proceed with measuring
#' @param max_strand_area, Maximum pixel area of a strand
#' @param cell_count Unique cell counter
#' @param KO_str string in filename corresponding to knockout genotype. Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype. Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout. Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout. Defaults to +/+.
#' @return Data frame with properties of synaptonemal (SC) measurements
#'
get_distance_general <- function(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,img_file,annotation, eccentricity_min, max_strand_area,cell_count,KO_str ,WT_str,KO_out, WT_out){
  tryCatch({
    no_strands <- nrow(num_strands)
    strand_count<- 0
    df_cols <- c("file","cell_id","genotype","strand_iter","foci_per_strand","iteration_on_strand","foci_x","foci_y","foci_x_line","foci_y_line","total_length", "distance_squared","distance_along","SC_pass_fail")
    dimensionless_dist <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
    while(strand_count<no_strands){
      strand_iter <- strand_iter + 1
      strand_count <- strand_count + 1
      # if area less than 150 pixels.. or not an outlier... keep
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
        ## once you have a single strand, get pixel info. Possibly save this? To delete later?
        single_info <- computeFeatures.shape(bwlabel(tmp_img))
        basic_info <- computeFeatures.basic(bwlabel(tmp_img),as.matrix(noise_gone))
        all_info <- computeFeatures(bwlabel(tmp_img),as.matrix(noise_gone))
        moment_info <- computeFeatures.moment(bwlabel(tmp_img),as.matrix(noise_gone))
        #moment_info <- computeFeatures.haralick(bwlabel(tmp_img),as.matrix(noise_gone))
        moment_info <- as.data.frame(moment_info)
        cx <- moment_info$m.cx
        cy <- moment_info$m.cy
        ## might actually want to find the real centre first..
        if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
          genotype <- WT_out
        }

        if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
          genotype <- KO_out
        }
        if (is.integer(nrow(per_strand_obj))){
          if(moment_info$m.eccentricity > eccentricity_min && nrow(per_strand_obj)>1){
            ## draw box around the middle
            ### find max, locally
            ### use noise_gone as original
            walkers <- 0*noise_gone
            ###
            noise_gone <- 2*noise_gone
            window <- 2
            window2 <- window*3
            ### start function here
            bright_loc <- find_start(window,noise_gone,cx,cy)
            mean_x <- as.numeric(bright_loc[1,1]) +cx -window2-1
            mean_y <- as.numeric(bright_loc[1,2]) +cy-window2-1
            ##
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
            ###########################################################################################
            ### first dir only deals with one half. But not a specific number...
            while(first_dir == 0){
              start_x <- round(mean(ix1))
              start_y <- round(mean(iy1))
              ### call get_first
              dir_1_out <- get_next_first_dir(new_square_1,ix1,iy1,dir_1,window,chosen_dir,distance_strand,first_dir)
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
              dir_2_out <- get_next_second_dir(new_square_2,ix2,iy2,dir_2,window,chosen_dir,distance_strand_2,second_dir)
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

            ### call measure distance between 2

            dimensionless_dist <- rbind(dimensionless_dist,get_distances_along(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_count,img_file,annotation,cell_count,strand_iter, per_strand_obj,KO_str ,WT_str,KO_out, WT_out))
            #dimensionless_dist <- get_distances_along(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_count,file,annotation,cell_count,strand_iter, per_strand_obj,KO_str ,WT_str,KO_out, WT_out)
            print("printing the 2 foci row")
            print(dimensionless_dist)

            ##### ends here
            ### you've got a single strand here. try and count distance between foci.
          }

          ### then the number of foci was 1
          else if(nrow(per_strand_obj)==1){
            dimensionless_dist <- rbind(dimensionless_dist,c(img_file,cell_count,genotype,strand_iter,1,1,"NA","NA","NA","NA","NA","NA", "NA","NA"))
            colnames(dimensionless_dist) <- df_cols
            print(dimensionless_dist)
            print("at strand number")
            print(strand_iter)

          }


        }
        ### then the number of foci was zero
        else{
          print("a strand with zero foci")
          dimensionless_dist <- rbind(dimensionless_dist,c(img_file,cell_count,genotype,strand_iter,0,"NA","NA","NA","NA","NA","NA","NA", "NA","NA"))
          colnames(dimensionless_dist) <- df_cols
          print(dimensionless_dist)
          print("at strand number")
          print(strand_iter)
        }

      }
    }
    return(dimensionless_dist)
  },

  error = function(e) {
    #what should be done in case of exception?
    str(e) # #prints structure of exception
    if(grepl( WT_str, img_file, fixed = TRUE) == TRUE){
      genotype <- WT_out
    }

    if(grepl( KO_str, img_file, fixed = TRUE) == TRUE){
      genotype <- KO_out
    }
    #dimensionless_dist_major_fail <- c(file, genotype, "NA", "NA", "NA", "fail")
    #return(dimensionless_dist_major_fail)

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
#' @param foci_label, black white mask with foci as objects, not necessarily on SC. Needed for computefeatures.
#' @param walkers, black white mask containing the line that traces through the middle of the SC. Computed earlier in get_distance
#' @param noise_gone, Black white mask with the single SC to be measured
#' @param start_x, x pixel location of where first branch terminated
#' @param start_y, y pixel location of where first branch terminated
#' @param start_x2, x pixel location of where second branch terminated
#' @param start_y2, y pixel location of where second branch terminated
#' @param start_dir, direction that the first branch was traveling along when it terminated
#' @param cx, centre of intensity x of the SC from computefeatures (annotation purposes only)
#' @param cy, centre of intensity y of the SC from computefeatures (annotation purposes only)
#' @param mean_x, starting point x that the two branches move away from to trace out the SC (annotation purposes only)
#' @param mean_y, starting point x that the two branches move away from to trace out the SC (annotation purposes only)
#' @param strand_iter, Strand number in iteration over all in cell
#' @param img_file, original filename that cell candidate came from. Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param cell_count Unique cell number
#' @param uid_strand Unique strand number
#' @param  per_strand_object Foci per strand object
#' @param KO_str string in filename corresponding to knockout genotype. Defaults to --.
#' @param WT_str string in filename corresponding to wildtype genotype. Defaults to ++.
#' @param KO_out string in output csv in genotype column, for knockout. Defaults to -/-.
#' @param WT_out string in output csv in genotype column, for knockout. Defaults to +/+.
#' @return List of fractional distances between foci for all SCs with two. Optional: total distances of SCs. Optional: images of all resulting traces/ foci locations.
#'
get_distances_along <- function(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_iter,img_file,annotation,cell_count, uid_strand,per_strand_object,KO_str ,WT_str,KO_out, WT_out){
  strand_info <- computeFeatures.moment(bwlabel(per_strand),as.matrix(foci_label))
  strand_info <- as.data.frame(strand_info)

  #### turn this into a table/ data frame
  no_foci <- nrow(strand_info)
  print("the number of foci in this strand is")
  print(no_foci)
  df_col <- c("file","cell_id","genotype","strand_iter","foci_per_strand","iteration_on_strand","foci_x","foci_y","foci_x_line","foci_y_line","total_length","distance_squared","distance_along","SC_pass_fail")
  foci_df <- data.frame(matrix(ncol = length(df_col), nrow = 0))
  my_walkers_matrix <- t(as.matrix(walkers))
  strand_test <- "pass"
  iter <- 0
  while(iter<= no_foci-1){
    iter <- iter + 1
    print("looking at position data of foci number")
    print(iter)
    foci_x <- strand_info$m.cx[iter]
    foci_y <- strand_info$m.cy[iter]
    #foci_df <- rbind(foci_df, c(foci_x,foci_y,0,0,0))
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
    #test_walker[x_curr,y_curr] <- 1

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
      ###
      colnames(foci_df) <- df_col
      #print("printing the x foci location")
      #print(foci_df$foci_x_line[iter] )
      #print("for iteration number")
      #print(iter)
      if(x_curr == as.numeric(foci_df$foci_x_line[iter]) && y_curr == as.numeric(foci_df$foci_y_line[iter])){
        ### record the distance along
        print("I found a focus, which is this many pixels along:")
        foci_dist_along <- dist_along
        print(dist_along)
        print("and its distance to the actual foci centre was")
        print(distance_fi)
        print("and this was a pass/fail")
        print(strand_test)
        if(annotation == "on"){
          ch1 <- bwlabel(noise_gone)
          ch1 <-channel(noise_gone,"grey")
          ch2 <- bwlabel(foci_label)
          ch2 <- channel(foci_label,"grey")
          bluered <- rgbImage(ch1, ch2, 0*ch1)
          plot(bluered)
          print("with foci locations included (actual)")
          plot(bluered)
          text(x = foci_x, y = foci_y, label = "+", col = "blue", cex = 2)
          print("with foci locations included, on the line")
          plot(bluered)
          text(x = mean_y, y = mean_x, label = "+", col = "magenta", cex = 2)
        }
        looping <- 0
        foci_df$distance_along[iter] <- dist_along
        foci_df$genotype[iter] <- genotype
        foci_df$SC_pass_fail[iter] <- strand_test

      }

      ### do stuff
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
        #}
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
    ### moving on to the next foci

  }
  colnames(foci_df) <- df_col

  if(length_walker<149){

    uid <- strand_iter
    # SIX new columns added: foci counter, foci per strand, f1 location x, f1 location y, f2 location x, f2 location y,  f1 location x (on line), f1 location y (on line), f2 location x (on line), f2 location y (on line),
    # df_cols <- c("file","genotype","foci_no","foci_per_strand", "total_SC_length","total_pixel_distance","foci_location_along", "fractional_distance_between_two", "pass_fail", "foci_location_x", "foci_location_y", "foci_location_x_line", "foci_location_y_line")
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
  else{
    #dimensionless_dist_fail_minor <- c(file, genotype, px_length,dim_length,(distance_strand+ distance_strand_2),"fail")
    #return(dimensionless_dist_fail_minor)
  }
  ## finish at start_x2

  ### checking the second arm
}


