#' measure_distances
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
#' @return Data frame with properties of synaptonemal (SC) measurements

# should take in same values as count_foci..
measure_distances <- function(img_path,offset_px = 0.2, offset_factor = 3, brush_size = 3, brush_sigma = 3, foci_norm = 0.01, annotation = "off", stage = "pachytene", eccentricity_min = 0.6, max_strand_area = 300)
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
  df_cols <- c("file","genotype","total_pixel_distance", "fractional_distance", "total_SC_length","pass_fail")
  df_lengths <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_lengths) <- df_cols
  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    filename_path_test = paste0(img_path,"/crops/",stage,"/", file)
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
      cell_count <- cell_count +1
      new_img<-img_orig
      #display(new_img)
      #### now see which have the right amount of strands
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)
      offset = offset_px
      thresh_crop = (new_img - localBackground > offset)
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
      foci_label <- threshold_foci_crop(img_orig_foci,offset_factor, brush_size, brush_sigma)
      ##### print properties of the images
      ### multiply strands by foci_label
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      coincident_foci <- bwlabel(foci_label*strands)
      overlap_no = table(coincident_foci)
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
      dimensionless_dist <- get_distance(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,file,annotation,eccentricity_min, max_strand_area)
      print("the dimension of the new row is")
      print(dim(dimensionless_dist))
      print("and part of the row is")
      print(dimensionless_dist)
      print("the full row is")
      df_lengths <- rbind(df_lengths, dimensionless_dist)
    }
  }
  colnames(df_lengths) <- df_cols
  return(df_lengths)

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
  disc = makeBrush(21, "disc")
  disc = disc / sum(disc)
  localBackground = filter2(image, disc)
  thresh_crop = (image - localBackground > offset)
  strands <- bwlabel(thresh_crop)
  return(strands)
}




#' threshold_foci_crop
#'
#' Creates mask for foci channel
#'
#' @param image foci channel image
#' @param offset_factor, Pixel value offset used in thresholding of foci channel
#' @param brush_size, size of brush to smooth the foci channel. Should be small to avoid erasing foci.
#' @param brush_sigma, sigma for Gaussian smooth of foci channel. Should be small to avoid erasing foci.
#' @return A black white mask with foci as objects
#'
threshold_foci_crop <- function(image, offset_factor, brush_size, brush_sigma){
  bg <- mean(image)
  offset = offset_factor*bg
  foci_th = image > bg + offset
  ### smooth it
  ### maybe up the contrast first??
  img_tmp_contrast = image
  #display(foci_mask_crop)
  w = makeBrush(size = brush_size, shape = 'gaussian', sigma = brush_sigma)
  img_flo = filter2(img_tmp_contrast, w)
  ## only choose objects above bright pixel value
  ## smooth foci channel
  foci_th = img_flo > bg + offset
  foci_label = bwlabel(foci_th)
  foci_label <- channel(foci_label, "grey")
  return(foci_label)
}



#' get_distance
#' Creates mask for SC channel
#'
#' @param strands, A black white mask with SCs as objects
#' @param num_strands, Number of individual strands on SC mask
#' @param new_img, Original strand/dna/SYCP3 channel image with noise removed.
#' @param foci_label, A black white mask with foci as objects
#' @param foci_count_strand, Number of foci counted located on the one SC
#' @param strand_iter, Strand number in iteration over all in cell
#' @param file, original filename that cell candidate came from. Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param eccentricity_min, The minimum eccentricity (from computefeatures) of a strand to proceed with measuring
#' @param max_strand_area, Maximum pixel area of a strand

#' @return Data frame with properties of synaptonemal (SC) measurements
#'
get_distance <- function(strands,num_strands,new_img,foci_label, foci_count_strand, strand_iter,file,annotation, eccentricity_min, max_strand_area){
  tryCatch({
    no_strands <- nrow(num_strands)
    strand_count<- 0
    while(strand_count<no_strands){
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
        if (is.integer(nrow(per_strand_obj))){
          if(moment_info$m.eccentricity > eccentricity_min && nrow(per_strand_obj) ==2){
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
            mean_x = as.numeric(bright_loc[1,1]) +cx -window2-1
            mean_y = as.numeric(bright_loc[1,2]) +cy-window2-1
            ##
            ix <- (round(mean_x)-window):(round(mean_x)+window)
            iy <- (round(mean_y)-window):(round(mean_y)+window)
            chosen_dir <- get_first_dir(noise_gone,ix,iy,window)
            walkers[round(mean(ix)),round(mean(iy))] = 1
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
            walkers[round(mean(ix1)),round(mean(iy1))] = 1
            walkers[round(mean(ix2)),round(mean(iy2))] = 1
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

              walkers[round(mean(ix1)),round(mean(iy1))] = 1
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
              distance_strand <- dir_1_out[(2*next_cord+2)]
              dir_1 <- dir_1_out[(2*next_cord+1)]
              second_dir <- dir_2_out[(2*next_cord+3)]
              start_dir2 <- dir_2_out[(2*next_cord+4)]
              walkers[round(mean(ix2)),round(mean(iy2))] = 1
              ## make the new cropped image.
              new_square_2 <-  noise_gone[ix2,iy2]
              if(distance_strand_2 >100){
                second_dir <- 1
              }
            }

            ### call measure distance between 2

            dimensionless_dist <- get_distance_between_two(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_iter,file,annotation)
            print(dimensionless_dist)
            return(dimensionless_dist)
            ##### ends here
            ### you've got a single strand here. try and count distance between foci.
          }
        }
        ## else: draw big square and find max.
        ###
      }
    }
  },

  error = function(e) {
    #what should be done in case of exception?
    str(e) # #prints structure of exception
    if(grepl( "++", file, fixed = TRUE) == TRUE){
      genotype <- "Fancm+/+"
    }

    if(grepl( "--", file, fixed = TRUE) == TRUE){
      genotype <- "Fancm-/-"
    }
    dimensionless_dist_major_fail <- c(file, genotype, "NA", "NA", "NA", "fail")
    return(dimensionless_dist_major_fail)

  }
  )

}



#' find_start
#'
#' Finds an appropriate place to start counting along an SC by looking in a window around the "centre of intensity"
#'
#' @param window window size that we want to look in
#' @param noise_gone Image containing only the SC of interest without background
#' @param cx,cy The location of the "centre of intensity" using computeFeatures (not necessarily on SC)
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
#' @param chosen_dir string: brightest direction of a line passing through starting point
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

#' @param new_square_1, the subsquare that contains the x y location of first branch walker in the middle.
#' @param ix1, current x position of first branch walker along SC
#' @param iy1, current y position of first branch walker along SC
#' @param dir_1, The direction (choice of 8) of the first branch step
#' @param window, number of pixels ahead that the intensity gradient is computed with
#' @param chosen_dir The brightest direction (choice of 4) of the previous step
#' @param distance_strand, current distance along the first branch
#' @param first_dir, zero while still measuring along first branch, one when first branch terminates and counting stops.
#' @return New sub square for first branch after taking one step
#'
get_next_first_dir <- function(new_square_1,ix1,iy1,dir_1,window,chosen_dir,distance_strand,first_dir){
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
  if(max_mean < 0.7){
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

  ## add to distance (1 or sqrt2 to distance calculation).

  ###########################################################################################
  ## if end of strand, set
  # eventually turn this into a function
  return(c(ix1,iy1,dir_1,distance_strand,first_dir,start_dir))

}





#' get_next_second_dir
#'
#' Moves one pixel away one second branch. Terminates if at the end of the SC.
#' @param new_square_2, the subsquare that contains the x y location of second branch walker in the middle.
#' @param ix2, current x position of second branch walker along SC
#' @param iy2, current y position of second branch walker along SC
#' @param dir_2, The direction (choice of 8) of the second (opposite to first) branch step
#' @param window, number of pixels ahead that the intensity gradient is computed with
#' @param chosen_dir The brightest direction (choice of 4) of the previous step
#' @param distance_strand_2, current distance along the second branch
#' @param second_dir, zero while still measuring along first branch, one when first branch terminates and counting stops.

#' @return New sub square for second branch after taking one step
#'

get_next_second_dir <- function(new_square_2,ix2,iy2,dir_2,window,chosen_dir,distance_strand_2,second_dir){
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
  ## need some condition to stop measuring? If mean is < 0.1? something like this..
  if(max_mean < 0.7){
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
#' get_distance_between_two
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
#' @param file, original filename that cell candidate came from. Used to identify e.g. genotype for data frame.
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @return List of fractional distances between foci for all SCs with two. Optional: total distances of SCs. Optional: images of all resulting traces/ foci locations.
#'
get_distance_between_two <- function(distance_strand,distance_strand_2,per_strand,foci_label, walkers, noise_gone,start_x,start_y,start_x2,start_y2,start_dir,cx,cy,mean_x,mean_y,strand_iter,file,annotation){
  print("we have a strand with two foci, located at")
  strand_info <- computeFeatures.moment(bwlabel(per_strand),as.matrix(foci_label))
  strand_info <- as.data.frame(strand_info)
  foci_1_x <- strand_info$m.cx[1]
  foci_1_y <- strand_info$m.cy[1]
  foci_2_x <- strand_info$m.cx[2]
  foci_2_y <- strand_info$m.cy[2]
  print(c(foci_1_x,foci_1_y,foci_2_x,foci_2_y))
  print("total distance is")
  print(distance_strand+ distance_strand_2)
  #### here is where you can identify the lengths.
  ### get the walkers matrix. Loop over, only if value = 1, assign a                    distance for a new matrix
  ### you've got a single strand here. try and count distance between foci.
  ### now loop over matrix
  my_walkers_matrix <- t(as.matrix(walkers))
  my_distance_matrix_f1 <- 0*my_walkers_matrix+100
  my_distance_matrix_f2 <- 0*my_walkers_matrix+100
  for(row in 1:nrow(my_walkers_matrix)) {
    for(col in 1:ncol(my_walkers_matrix)) {
      if(my_walkers_matrix[row, col]==1){
        my_distance_matrix_f1[row, col] <- (row-foci_1_y)^2+(col-foci_1_x)^2
        my_distance_matrix_f2[row, col] <- (row-foci_2_y)^2+(col-foci_2_x)^2
      }
    }
  }


  #### find max. plot these onto the images.
  bright_loc_f1 <- which(my_distance_matrix_f1 == min(my_distance_matrix_f1),arr.ind = TRUE)
  mean_x_f1 = as.numeric(bright_loc_f1[1,1])
  mean_y_f1 = as.numeric(bright_loc_f1[1,2])

  distance_f1 <- (foci_1_y-mean_x_f1)^2 +(foci_1_x-mean_y_f1)^2

  ###
  bright_loc_f2 <- which(my_distance_matrix_f2 == min(my_distance_matrix_f2),arr.ind = TRUE)
  mean_x_f2 = as.numeric(bright_loc_f2[1,1])
  mean_y_f2 = as.numeric(bright_loc_f2[1,2])


  distance_f2 <- (foci_2_y-mean_x_f2)^2 +(foci_2_x-mean_y_f2)^2
  # deleting for now
  if (annotation=="on"){
    ch1 = bwlabel(walkers)
    ch1 <- channel(ch1, "grey")
    ch2 = bwlabel(noise_gone)
    ch2 <-channel(noise_gone,"grey")
    ch3 = bwlabel(per_strand)
    ch3 <- channel(per_strand,"grey")
    bluered <- rgbImage(ch2, ch1, ch3)
    plot(bluered)
    text(x = foci_1_x, y = foci_1_y, label = "+", col = "yellow", cex = 2)
    text(x = foci_2_x, y = foci_2_y, label = "+", col = "yellow", cex = 2)
    text(x = mean_y_f1, y = mean_x_f1, label = "+", col = "magenta", cex = 2)
    text(x = mean_y_f2, y = mean_x_f2, label = "+", col = "green", cex = 2)
    text(x = start_x2, y = start_y2, label = "+", col = "blue", cex = 2)
    text(x = start_x, y = start_y, label = "+", col = "blue", cex = 2)

  }

  ### now that you have the exact positions of the foci on the walker channel, we want to count distances.

  ### start at start_x
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
  test_walker[x_curr,y_curr] <- 1

  ## measuring foci
  dist_between_foci <- 0
  measuring_distance <- 0
  foci_out_2 <- 0
  #while(x_curr != start_x2 & y_curr != start_y2){
  #while(length_walker < 100){
  while(looping == 1){
    length_walker <- length_walker + 1


    if(x_curr == start_x2){
      if(y_curr == start_y2){
        looping <- 0
      }
    }

    if (length_walker > 150){
      print("probably an infinite loop")
      looping <- 0
    }



    ###
    if(x_curr == mean_y_f1 && y_curr == mean_x_f1 | x_curr == mean_y_f2 && y_curr == mean_x_f2 ){
      foci_out_2 <- foci_out_2+1
    }

    if(foci_out_2 == 1){
      measuring_distance <- 1
    }

    if(foci_out_2 > 1){
      measuring_distance <- 0
    }

    test_walker[x_curr,y_curr] <- 1
    ### do stuff
    ### do the opposite of the travelling direction

    if(dir_curr == 1){
      # 8 position
      if(walkers[x_curr-1,y_curr-1]==1){
        y_curr <- y_curr -1
        x_curr <- x_curr -1
        dir_curr <- 8
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 1 position
      else if(walkers[x_curr,y_curr-1]==1){
        y_curr <- y_curr -1
        dir_curr <- 1
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }

      # 2 position
      else if(walkers[x_curr+1,y_curr-1]==1){
        y_curr <- y_curr -1
        x_curr <- x_curr +1
        dir_curr <- 2
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }

    }

    else if(dir_curr == 2){
      #}
      # 1 position
      if(walkers[x_curr,y_curr-1]==1){
        y_curr <- y_curr -1
        dir_curr <- 1
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }

      # 2 position
      else if(walkers[x_curr+1,y_curr-1]==1){
        y_curr <- y_curr -1
        x_curr <- x_curr +1
        dir_curr <- 2
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }

      }

      # 3 position
      else if(walkers[x_curr+1,y_curr]==1){
        x_curr <- x_curr +1
        dir_curr <- 3
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }

    }
    else if(dir_curr == 3){
      # 2 position
      if(walkers[x_curr+1,y_curr-1]==1){
        y_curr <- y_curr -1
        x_curr <- x_curr +1
        dir_curr <- 2
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }

      # 3 position
      else if(walkers[x_curr+1,y_curr]==1){
        x_curr <- x_curr +1
        dir_curr <- 3
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 4 position
      else if(walkers[x_curr+1,y_curr+1]==1){
        y_curr <- y_curr +1
        x_curr <- x_curr +1
        dir_curr <- 4
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }

    }
    else if(dir_curr == 4){
      # 3 position
      if(walkers[x_curr+1,y_curr]==1){
        x_curr <- x_curr +1
        dir_curr <- 3
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 4 position
      else if(walkers[x_curr+1,y_curr+1]==1){
        y_curr <- y_curr +1
        x_curr <- x_curr +1
        dir_curr <- 4
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 5 position
      else if(walkers[x_curr,y_curr+1]==1){
        y_curr <- y_curr +1
        dir_curr <- 5
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
    }
    else if(dir_curr == 5){
      # 4 position
      if(walkers[x_curr+1,y_curr+1]==1){
        y_curr <- y_curr +1
        x_curr <- x_curr +1
        dir_curr <- 4
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 5 position
      else if(walkers[x_curr,y_curr+1]==1){
        y_curr <- y_curr +1
        dir_curr <- 5
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 6 position
      else if(walkers[x_curr-1,y_curr+1]==1){
        x_curr <- x_curr -1
        y_curr <- y_curr +1
        dir_curr <- 6
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
    }
    else if(dir_curr == 6){
      # 5 position
      if(walkers[x_curr,y_curr+1]==1){
        y_curr <- y_curr +1
        dir_curr <- 5
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 6 position
      else if(walkers[x_curr-1,y_curr+1]==1){
        x_curr <- x_curr -1
        y_curr <- y_curr +1
        dir_curr <- 6
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 7 position
      else if(walkers[x_curr-1,y_curr]==1){
        x_curr <- x_curr -1
        dir_curr <- 7
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }

    }
    else if(dir_curr == 7){
      # 6 position
      if(walkers[x_curr-1,y_curr+1]==1){
        x_curr <- x_curr -1
        y_curr <- y_curr +1
        dir_curr <- 6
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 7 position
      else if(walkers[x_curr-1,y_curr]==1){
        x_curr <- x_curr -1
        dir_curr <- 7
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 8 position
      else if(walkers[x_curr-1,y_curr-1]==1){
        x_curr <- x_curr -1
        y_curr <- y_curr -1
        dir_curr <- 8
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }

    }
    else if(dir_curr == 8){
      # 7 position
      if(walkers[x_curr-1,y_curr]==1){
        x_curr <- x_curr -1
        dir_curr <- 7
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
      # 8 position
      else if(walkers[x_curr-1,y_curr-1]==1){
        x_curr <- x_curr -1
        y_curr <- y_curr -1
        dir_curr <- 8
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + sqrt(2)
        }
      }
      # 1 position
      else if(walkers[x_curr,y_curr-1]==1){
        y_curr <- y_curr -1
        dir_curr <- 1
        length_walker <- length_walker + 1
        if(measuring_distance == 1){
          dist_between_foci <- dist_between_foci + 1
        }
      }
    }


    ### find direction to move i.e. find the closest 1.
    ## just look at the whole crop....
    ### stop doing stuff
  }
  if(annotation == "on"){
    print("here are the results for this SC")
    display(rgbImage(walkers, test_walker, test_walker))
    # deleting for now
    text(x = mean_y_f1, y = mean_x_f1, label = "+", col = "magenta", cex = 2)
    text(x = mean_y_f2, y = mean_x_f2, label = "+", col = "magenta", cex = 2)

  }


  dim_length <- dist_between_foci/(distance_strand+ distance_strand_2)
  px_length <- dist_between_foci

  if(length_walker<149){
    strand_iter <- strand_iter +1
    print("on iteration")
    print(strand_iter)
    if (dim_length >1e-6 && dim_length < 1 && (distance_strand+ distance_strand_2) > 0){
      if (foci_out_2 >1){
        if(distance_f1 < 10){
          if(distance_f2 < 10){
            #dimensionless_dist <- append(dimensionless_dist,px_length)
            if(grepl( "++", file, fixed = TRUE) == TRUE){
              genotype <- "Fancm+/+"
            }

            if(grepl( "--", file, fixed = TRUE) == TRUE){
              genotype <- "Fancm-/-"
            }
            dimensionless_dist_pass <- c(file, genotype, px_length,dim_length,(distance_strand+ distance_strand_2),"pass")
            print("This strand managed to pass through:")
            ch1 = bwlabel(walkers)
            ch1 <- channel(ch1, "grey")
            ch2 = bwlabel(noise_gone)
            ch2 <-channel(noise_gone,"grey")
            ch3 = bwlabel(per_strand)
            ch3 <- channel(per_strand,"grey")
            bluered <- rgbImage(ch2, ch1, ch3)
            #print("break")
            #display(bluered)
            #bluered <- rgbImage(ch2, ch1, ch1)
            plot(bluered)
            text(x = foci_1_x, y = foci_1_y, label = "+", col = "yellow", cex = 2)
            text(x = foci_2_x, y = foci_2_y, label = "+", col = "yellow", cex = 2)
            text(x = mean_y_f1, y = mean_x_f1, label = "+", col = "magenta", cex = 2)
            text(x = mean_y_f2, y = mean_x_f2, label = "+", col = "green", cex = 2)
            text(x = start_x2, y = start_y2, label = "+", col = "blue", cex = 2)
            text(x = start_x, y = start_y, label = "+", col = "blue", cex = 2)
            print("This one worked. check that it walked successfully")
            display(rgbImage(walkers, test_walker, test_walker))
            # deleting for now
            text(x = mean_y_f1, y = mean_x_f1, label = "+", col = "magenta", cex = 2)
            text(x = mean_y_f2, y = mean_x_f2, label = "+", col = "magenta", cex = 2)
            return(dimensionless_dist_pass)

          }
        }

      }


    }


  }
  else{
    if(annotation == "on"){
      print("the following failed and will be excluded")

      plot(noise_gone)
      text(x = foci_1_x, y = foci_1_y, label = "+", col = "red", cex = 2)
      text(x = foci_2_x, y = foci_2_y, label = "+", col = "blue", cex = 2)
      ch1 = bwlabel(walkers)
      ch1 <- channel(ch1, "grey")
      ch2 = bwlabel(noise_gone)
      ch2 <-channel(noise_gone,"grey")
      ch3 = bwlabel(per_strand)
      ch3 <- channel(per_strand,"grey")
      bluered <- rgbImage(ch2, ch1, ch3)
      #print("break")
      #display(bluered)
      bluered <- rgbImage(ch2, ch1, ch1)
      #print("break")
      #display(bluered)
      print(file)

      ### deleting for now
      bluered <- rgbImage(ch2, ch3, 0*ch3)
      plot(bluered)
      text(x = foci_1_x, y = foci_1_y, label = "+", col = "yellow", cex = 2)
      text(x = foci_2_x, y = foci_2_y, label = "+", col = "yellow", cex = 2)
      text(x = cx, y = cy, label = "+", col = "blue", cex = 2)
      text(x = mean_x, y = mean_y, label = "+", col = "magenta", cex = 2)



    }
    if(grepl( "++", file, fixed = TRUE) == TRUE){
      genotype <- "Fancm+/+"
    }

    if(grepl( "--", file, fixed = TRUE) == TRUE){
      genotype <- "Fancm-/-"
    }
    dimensionless_dist_fail_minor <- c(file, genotype, px_length,dim_length,(distance_strand+ distance_strand_2),"fail")
    return(dimensionless_dist_fail_minor)
  }
  ## finish at start_x2

  ### checking the second arm
}


