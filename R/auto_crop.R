#' auto_crop
#'
#' crop an image around each viable cell candidate.
#' @importFrom stats median sd
#' @importFrom EBImage bwlabel channel colorLabels computeFeatures computeFeatures.basic computeFeatures.moment computeFeatures.shape display filter2 makeBrush readImage rgbImage rmObjects rotate writeImage
#' @importFrom graphics text
#' @importFrom utils str
#' @export auto_crop
#' @param img_path, path containing image data to analyse
#' @param max_cell_area, Maximum pixel area of a cell candidate
#' @param min_cell_area, Minimum pixel area of a cell candidate
#' @param mean_pix, Mean pixel intensity of cell crop (in SYCP3 channel) for normalisation
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param blob_factor, Contrast factor to multiply original image by before smoothing/smudging
#' @param bg_blob_factor, Contrast factor to multiply original image by to take background. Used prior to thresholding.
#' @param offset, Pixel value offset from bg_blob_factor. Used in thresholding to make blob mask.
#' @param final_blob_amp, Contrast factor to multiply smoothed/smudged image. Used in thresholding to make blob mask.
#' @param test_amount, Optional number of first N images you want to run function on. For troubleshooting/testing/variable calibration purposes.
#' @param brush_size_blob, Brush size for smudging the dna channel to make blobs
#' @param sigma_blob, Sigma in Gaussian brush for smudging the dna channel to make blobs
#' @param cell_aspect_ratio Maximum aspect ratio of blob to be defined as a cell

#' @return cropped SC and foci channels around single cells, regardless of stage


auto_crop <- function(img_path,  max_cell_area = 20000, min_cell_area = 7000, mean_pix = 0.08, annotation = "off", blob_factor = 15, bg_blob_factor = 10,  offset = 0.2, final_blob_amp = 10, test_amount = 0,brush_size_blob = 51, sigma_blob = 15, cell_aspect_ratio = 2)
{
  file_list <- list.files(img_path)
  dir.create(paste0(img_path,"/crops"))
  # input :

  # output : a bunch of output jpegs? Or save them all?

  #BiocManager::install("EBImage")
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  antibody3_store <- 0

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    print(file)
    file_base = file
    filename_path_test = paste0(img_path,"/", file)
    file = filename_path_test
    if(grepl("*DAPI.jpeg$", file)){
      file_DAPI = file
      image <- readImage(file_DAPI)
      img_orig_DAPI <- channel(image, "grey")
      antibody3_store <- 1
    }
    if(grepl("*SYCP3.jpeg$", file)){

      file_dna = file

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
      image_count <- image_count +1
      if(test_amount != 0 && image_count > test_amount){
        break

      }

      #### function: blur the image
      ## call it on img_orig, optional offset
      blob_th <- get_blobs(img_orig,blob_factor, bg_blob_factor, offset,final_blob_amp,brush_size_blob, sigma_blob)

      blob_label = bwlabel(blob_th)
      blob_label <- channel(blob_label, "gray")
      candidate <- bwlabel(blob_label)

      ## function: remove things that aren't cells

      retained <- keep_cells(candidate, max_cell_area, min_cell_area,cell_aspect_ratio)

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
        crop_single_object(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_dna,file_foci,file_DAPI,cell_count, mean_pix, annotation, file_base, img_path)

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
#' @export
#' @param img_orig Original image
#'
#' @param blob_factor, Contrast factor to multiply original image by before smoothing/smudging
#' @param bg_blob_factor, Contrast factor to multiply original image by to take background. Used prior to thresholding.
#' @param offset, Pixel value offset from bg_blob_factor. Used in thresholding to make blob mask.
#' @param final_blob_amp, Contrast factor to multiply smoothed/smudged image. Used in thresholding to make blob mask.
#' @param brush_size_blob, Brush size for smudging the dna channel to make blobs
#' @param sigma_blob, Sigma in Gaussian brush for smudging the dna channel to make blobs
#' @return Mask with cell candidates

#################################### new function ####################################
get_blobs <- function(img_orig, blob_factor, bg_blob_factor, offset,final_blob_amp, brush_size_blob,sigma_blob){

  # input:

  # output:

  ## subfunction: get signals and make BW mask
  ### default offset
  thresh <- blob_factor*img_orig
  # subfunction: big blur to blobs
  img_tmp_dna <- img_orig
  img_tmp <- thresh
  #w = makeBrush(size = 51, shape = 'gaussian', sigma = 15)
  w = makeBrush(size = brush_size_blob, shape = 'gaussian', sigma = sigma_blob)
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
#' @export
#' @param candidate Mask of individual cell candidates
#' @param max_cell_area, Maximum pixel area of a cell candidate
#' @param min_cell_area, Minimum pixel area of a cell candidate
#' @param cell_aspect_ratio Maximum aspect ratio of blob to be defined as a cell

#' @return Mask of cell candidates which meet size criteria
keep_cells <- function(candidate, max_cell_area, min_cell_area, cell_aspect_ratio){

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
    if(semi_maj/semi_min > cell_aspect_ratio & is.na(semi_maj/semi_min)==FALSE){
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
#' @export
#' @param retained Mask of cell candidates which meet size criteria. After smoothing/smudging and thresholding.
#' @param OOI_final, Objects of interest count. Total number of cell candidates in retained.
#' @param counter_final, Counter for single cell we are focussing on. Remove all other cells where counter_single not equal to counter_final.
#' @param img_orig, description
#' @param img_orig_foci, description
#' @param img_orig_DAPI, description
#' @param file_dna, filename of dna channel image
#' @param file_foci, filename of foci channel image
#' @param file_DAPI, filename of DAPI channel image
#' @param cell_count, counter for successful crops around cells
#' @param mean_pix, Mean pixel intensity of cell crop (in SYCP3 channel) for normalisation
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param file_base, filename base common to all three channels i.e. without -MLH3.jpeg etc.
#' @param img_path, path containing image data to analyse

#'



#' @return Crops around all candidates in both channels
#'
crop_single_object <- function(retained, OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_dna,file_foci,file_DAPI,cell_count, mean_pix, annotation, file_base, img_path){
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
  my_matrix = t(as.matrix(xx))
  i <- 0
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
      print(file_dna)
      file_dna <- gsub('-SYCP3.jpeg','', file_base)
      print(file_dna)
      filename_crop = paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,"-SYCP3.jpeg")
      writeImage(new_img, filename_crop)

      new_img_foci <- noise_gone_foci[ix, iy]
      #file_foci <- tools::file_path_sans_ext(file_foci)
      file_foci <- gsub('-MLH3.jpeg','', file_foci)
      filename_crop_foci = paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,"-MLH3.jpeg")
      writeImage(new_img_foci, filename_crop_foci)

      new_img_DAPI <- noise_gone_DAPI[ix, iy]
      #file_foci <- tools::file_path_sans_ext(file_foci)
      file_DAPI <- gsub('-DAPI.jpeg','', file_DAPI)
      filename_crop_DAPI = paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,"-DAPI.jpeg")
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

