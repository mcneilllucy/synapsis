#' auto_crop
#'
#' crop an image around each viable cell candidate.
#'
#' @import EBImage
#' @import stats
#' @import graphics
#' @import utils
#' @export auto_crop
#' @param file_list The file list
#' @param img_path The path
#' @return cropped SC and foci channels around single cells, regardless of stage


auto_crop <- function(img_path,  max_cell_area = 20000, min_cell_area = 7000, mean_pix = 0.08, annotation = "off", blob_factor = 15, bg_blob_factor = 10,  offset = 0.2, final_blob_amp = 10)
{
  file_list <- list.files(img_path)
  setwd(img_path)
  dir.create("crops")
  # input :

  # output : a bunch of output jpegs? Or save them all?

  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  antibody3_store <- 0

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    setwd(img_path)
    if(grepl("*DAPI.jpeg$", file)){
      file_DAPI = file
      image <- readImage(file_DAPI)
      img_orig_DAPI <- channel(image, "grey")
      antibody3_store <- 1
    }
    if(grepl("*SYCP3.jpeg$", file)){

      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(2*image, "grey")
      antibody1_store <- 1
    }
    if(grepl("*MLH3.jpeg$", file)){
      file_foci = file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(antibody1_store + antibody2_store + antibody3_store ==3){

      #### function: blur the image
      ## call it on img_orig, optional offset
      blob_th <- get_blobs(img_orig,blob_factor, bg_blob_factor, offset,final_blob_amp)

      blob_label = bwlabel(blob_th)
      blob_label <- channel(blob_label, "gray")
      candidate <- bwlabel(blob_label)

      ## function: remove things that aren't cells

      retained <- keep_cells(candidate, max_cell_area, min_cell_area)

      ### crop over each cell
      ## function: crop around every single viable cell
      ### crop foci channel here
      ## Loop over all objects in this final "removed" image


      x_final <- computeFeatures.shape(retained)
      x_final <- data.frame(x_final)
      OOI_final <- nrow(x_final)

      counter_final <- 0
      # looping through each object to crop
      while(counter_final<OOI_final){
        counter_final <- counter_final+1


        cell_count <- cell_count +1
        crop_single_object(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_dna,file_foci,file_DAPI,cell_count, mean_pix, annotation)

        print("Crop number:")
        print(cell_count)
      }
      antibody1_store <- 0
      antibody2_store <- 0
      antibody3_store <- 0
    }

  }

print("out of")
print(image_count)
print("images, we got")
print(cell_count)
print("viable cells")
}


### add all the functions used here

#' get_blobs
#'
#' Makes mask of all objects bright enough
#'
#' @import EBImage
#' @export
#' @param img_orig Original image
#' @return Mask with cell candidates

#################################### new function ####################################
get_blobs <- function(img_orig, blob_factor, bg_blob_factor, offset,final_blob_amp){

  # input:

  # output:

  ## subfunction: get signals and make BW mask
  ### default offset
  thresh <- blob_factor*img_orig
  # subfunction: big blur to blobs
  img_tmp_dna <- img_orig
  img_tmp <- thresh
  w = makeBrush(size = 51, shape = 'gaussian', sigma = 15)
  img_flo = filter2(img_tmp, w)
  ## default amplification
  bg <- mean(bg_blob_factor*img_tmp)
  #if(crop_method == "regular"){
  blob_th = final_blob_amp*img_flo > bg + offset
  #}

  #if(crop_method == "watershed"){
  #  blob_th = 10*img_flo > bg + offset
  #}
  return(blob_th)

}

#################################### new function ####################################


#' keep_cells
#'
#' Deletes objects in mask which are too small, large, oblong i.e. unlikely to be a cell
#'
#' @import EBImage
#' @export
#' @param candidate Mask of individual cell candidates
#' @return Mask of cell candidates which meet size criteria
keep_cells <- function(candidate, max_cell_area, min_cell_area){


  # input:

  # output:



  # delete everything that's too small
  colorimg<- colorLabels(candidate, normalize = TRUE)
  x <- computeFeatures.shape(candidate)
  x <- data.frame(x)
  OOI <- nrow(x)
  counter <- 0
  retained <- candidate

  while(counter<OOI){
    counter <- counter+1
    pixel_area = x$s.area[counter]
    semi_maj <- x$s.radius.max[counter]
    semi_min <- x$s.radius.min[counter]
    # if statement checking if it's the wrong area
    if(pixel_area> max_cell_area | pixel_area < min_cell_area){
      retained <- as.numeric(retained)*rmObjects(candidate, counter, reenumerate = TRUE)
    }
    ## if statement checking that it's not too long i.e. not at edge.
    if(semi_maj/semi_min > 2 & is.na(semi_maj/semi_min)==FALSE){
      retained <- as.numeric(retained)*rmObjects(candidate, counter, reenumerate = TRUE)
    }

  }
  retained <- bwlabel(retained)
  return(retained)
}

#################################### new function ####################################


#' crop_single_object
#'
#' Creates mask for every individual cell candidate in mask
#'
#' @import EBImage
#' @export
#' @param retained Mask of cell candidates which meet size criteria
#' @return Crops aroudn all candidates in both channels
#'
crop_single_object <- function(retained, OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_dna,file_foci,file_DAPI,cell_count, mean_pix, annotation){
  tmp_img <- retained
  ## have a single object
  ### delete all other objects
  counter_single <- 0
  # looping over all other objects to crop
  while(counter_single < OOI_final){
    counter_single <- counter_single + 1
    # iteratively remove all other objects
    if(counter_single != counter_final){
      tmp_img <- as.numeric(tmp_img)*rmObjects(bwlabel(retained), counter_single, reenumerate = TRUE)

    }


  }
  if(annotation == "on"){
    print("here is the mask of a single cell")
    display(tmp_img)

  }

  ## function: remove noise
  noise_gone <- bwlabel(tmp_img)*as.matrix(img_orig)
  noise_gone_foci <- bwlabel(tmp_img)*as.matrix(img_orig_foci)
  noise_gone_DAPI <- bwlabel(tmp_img)*as.matrix(img_orig_DAPI)
  ## first get the row and column list that has a one in it.
  row_list <- c()
  col_list <- c()
  # I think this is quick enough for now.. takes less than 10s...
  xx = data.frame(as.numeric(tmp_img))
  xx <- data.frame(bwlabel(tmp_img))
  xm2 = t(as.matrix(xx))
  i <- 0
  my_matrix <- xm2
  ### now loop over matrix
  for(row in 1:nrow(my_matrix)) {
    for(col in 1:ncol(my_matrix)) {
      if(my_matrix[row, col]==1){
        row_list[i] <- row
        col_list[i] <- col
        i <- i+1
      }
    }
  }

  if (length(col_list>1)){
    cy <- mean(row_list)
    cx <- mean(col_list)

    x_maj <- max(col_list)
    x_min <- min(col_list)

    y_maj <- max(row_list)
    y_min <- min(row_list)
    ## radius
    ## total number of pixels given by list length. Find radius of area.
    ## crops to a square for the moment. can change this.
    max_r <- max(x_maj-cx,y_maj-cy)
    crop_r <- floor(max_r)
    # might want to do this as a matrix
    top_left_x <- floor(cx-crop_r)
    top_left_y <- floor(cy-crop_r)

    bottom_left_x <- floor(cx-crop_r)
    bottom_left_y <- floor(cy+crop_r)

    bottom_right_x <- floor(cx+crop_r)
    bottom_right_y <- floor(cy+crop_r)

    top_right_x <- floor(cx-crop_r)
    top_right_y <- floor(cy+crop_r)

    ## crop image
    ix <- bottom_left_x:bottom_right_x
    iy <- top_left_y:bottom_left_y

    # determine the dimensions, 2 in this case

    ## cropping part
    tryCatch({
      new_img <- noise_gone[ix, iy]
      ## want all images to have the same mean to 0.1
      orig_mean <- mean(new_img)
      mean_factor <- mean_pix/orig_mean
      new_img <- new_img*mean_factor
      #file_dna <- tools::file_path_sans_ext(file_dna)
      file_dna <- gsub('-SYCP3.jpeg','', file_dna)
      filename_crop = paste0("./crops/", file_dna,"-crop-",cell_count,"-SYCP3.jpeg")
      writeImage(new_img, filename_crop)

      new_img_foci <- noise_gone_foci[ix, iy]
      #file_foci <- tools::file_path_sans_ext(file_foci)
      file_foci <- gsub('-MLH3.jpeg','', file_foci)
      filename_crop_foci = paste0("./crops/", file_foci,"-crop-",cell_count,"-MLH3.jpeg")
      writeImage(new_img_foci, filename_crop_foci)

      new_img_DAPI <- noise_gone_DAPI[ix, iy]
      #file_foci <- tools::file_path_sans_ext(file_foci)
      file_DAPI <- gsub('-DAPI.jpeg','', file_DAPI)
      filename_crop_DAPI = paste0("./crops/", file_DAPI,"-crop-",cell_count,"-DAPI.jpeg")
      writeImage(new_img_DAPI, filename_crop_DAPI)

      if(annotation=="on"){
        print("from the file:")
        print(file_dna)
        display(img_orig)
        print("I cropped this cell:")
        display(new_img)
        print("whose cell number is")
        print(cell_count)
      }


      #### strand related stuff here
    },
    error = function(e) {
      #what should be done in case of exception?
      str(e) # #prints structure of exception
      print("couldn't crop it")

    }
    )
  }

}


#################################### new function ####################################

