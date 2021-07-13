#' auto_crop_fast
#'
#' crop an image around each viable cell candidate.
#' @importFrom stats median sd
#' @importFrom EBImage bwlabel channel colorLabels computeFeatures computeFeatures.basic computeFeatures.moment computeFeatures.shape display filter2 makeBrush readImage rgbImage rmObjects rotate writeImage watershed
#' @importFrom graphics text
#' @importFrom utils str
#' @export auto_crop_fast
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
#' @param cell_aspect_ratio Maximum aspect ratio of blob to be defined as a cell
#' @param sigma_blob, Sigma in Gaussian brush for smudging the dna channel to make blobs
#' @param channel1_string String appended to the files showing the channel illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel illuminating synaptonemal complexes. Defaults to SYCP3
#' @param channel3_string Optional. String appended to the files showing the channel illuminating cell structures. Defaults to DAPI, if third channel == "on".
#' @param third_channel Optional, defaults to "off".
#' @param file_ext file extension of your images e.g. tif jpeg or png.
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' auto_crop_fast(demo_path, annotation = "on", max_cell_area = 30000,
#' min_cell_area = 7000)
#' @author Lucy McNeill
#' @return cropped SC and foci channels around single cells, regardless of stage

auto_crop_fast <- function(img_path,  max_cell_area = 20000, min_cell_area = 7000, mean_pix = 0.08, annotation = "off", blob_factor = 15, bg_blob_factor = 10,  offset = 0.2, final_blob_amp = 10, test_amount = 0,brush_size_blob = 51, sigma_blob = 15, channel3_string = "DAPI", channel2_string = "SYCP3", channel1_string = "MLH3", file_ext = "jpeg", third_channel = "off",cell_aspect_ratio = 2)
{
  file_list <- list.files(img_path)
  dir.create(paste0(img_path,"/crops"))
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  antibody3_store <- 0
  ## for each image that is *-dna.jpeg,
  for (img_file in file_list){
    file_base <- img_file
    filename_path_test <- paste0(img_path,"/",img_file)
    img_file <- filename_path_test
    if(third_channel == "on"){
      if(grepl(paste0('*',channel3_string,'.',file_ext,'$'),img_file)){
        file_DAPI <- img_file
        image <- readImage(file_DAPI)
        img_orig_DAPI <- channel(image, "grey")
        antibody3_store <- 1
      }
    }
    if(grepl(paste0('*',channel2_string,'.',file_ext,'$'), img_file)){
      file_dna <- img_file
      print(file)
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
    if(third_channel == "off"){
      antibody3_store <- 1
    }
    if(antibody1_store + antibody2_store + antibody3_store ==3){
      image_count <- image_count +1
      if(test_amount != 0 && image_count > test_amount){
        break
      }
      #### blur the image with get_blobs. call it on img_orig, optional offset
      blob_th <- get_blobs(img_orig,blob_factor, bg_blob_factor, offset,final_blob_amp,brush_size_blob, sigma_blob)
      blob_label <- bwlabel(blob_th)
      blob_label <- channel(blob_label, "gray")
      candidate <- bwlabel(blob_label)
      ## remove things that aren't cells
      retained <- keep_cells(candidate, max_cell_area, min_cell_area,cell_aspect_ratio)
      ### crop over each cell
      ## Loop over all objects in this final "removed" image
      x_final <- computeFeatures.shape(retained)
      moment_info <- computeFeatures.moment(retained,as.matrix(img_orig))
      x_final <- data.frame(x_final)
      moment_info <- as.data.frame(moment_info)
      OOI_final <- nrow(x_final)
      counter_final <- 0
      r_max <- x_final$s.radius.max
      cx <- moment_info$m.cx
      cy <- moment_info$m.cy
      # looping through each object to crop
      while(counter_final<OOI_final){
        counter_final <- counter_final+1
        ### row of interest is the counter_final'th row of x_final
        cell_count <- cell_count +1
        if(third_channel=="on"){
          crop_single_object_fast(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_dna,file_foci,file_DAPI,cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel)
        }
        else{
          crop_single_object_fast(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_foci,file_dna,file_foci,file_foci,cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel)
        }
      }
      antibody1_store <- 0
      antibody2_store <- 0
      antibody3_store <- 0
    }
  }
print("out of")
print(image_count-1)
print("images, we got")
print(cell_count)
print("viable cells")
}


#' crop_single_object_fast
#'
#' Creates mask for every individual cell candidate in mask
#'
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
#' @param r_max maximum radius of blob for cropping
#' @param cx centre of blob x
#' @param cy centre of blob y
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param file_base, filename base common to all three channels i.e. without -MLH3.jpeg etc.
#' @param img_path, path containing image data to analyse
#' @param channel1_string String appended to the files showing the channel illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel illuminating synaptonemal complexes. Defaults to SYCP3
#' @param channel3_string Optional. String appended to the files showing the channel illuminating cell structures. Defaults to DAPI, if third channel == "on".
#' @param third_channel Optional, defaults to "off".
#' @param file_ext file extension of your images e.g. tif jpeg or png.



#' @return Crops around all candidates in both channels
#'
crop_single_object_fast <- function(retained, OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI="blank",file_dna,file_foci,file_DAPI = "blank",cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel){
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
  noise_gone <- bwlabel(tmp_img)*as.matrix(img_orig)
  noise_gone_foci <- bwlabel(tmp_img)*as.matrix(img_orig_foci)
  if(third_channel == "on"){
    noise_gone_DAPI <- bwlabel(tmp_img)*as.matrix(img_orig_DAPI)
  }
  ### use the features of tmp_img
  ### here we have the single object. Need to identify its centre value and radius..
  crop_r <- floor(r_max[counter_final])
  cx <- cx[counter_final]
  cy <- cy[counter_final]
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
### cropping finished
  ## cropping part
  tryCatch({
    new_img <- noise_gone[ix, iy]
    ## want all images to have the same mean (mean_pix)
    orig_mean <- mean(new_img)
    mean_factor <- mean_pix/orig_mean
    new_img <- new_img*mean_factor
    file_stub <- paste0('-',channel2_string,'.',file_ext)
    file_dna <- gsub(file_stub,'', file_base)
    filename_crop <- paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,file_stub)
    writeImage(new_img, filename_crop)
    new_img_foci <- noise_gone_foci[ix, iy]
    file_stub <- paste0('-',channel1_string,'.',file_ext)
    file_foci <- gsub(file_stub,'', file_foci)
    filename_crop_foci <- paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,file_stub)
    writeImage(new_img_foci, filename_crop_foci)
    if(third_channel == "on"){
      new_img_DAPI <- noise_gone_DAPI[ix, iy]
      file_stub <- paste0('-',channel3_string,'.',file_ext)
      file_DAPI <- gsub(file_stub,'', file_DAPI)
      filename_crop_DAPI <- paste0(img_path,"/crops/", file_dna,"-crop-",cell_count,file_stub)
      writeImage(new_img_DAPI, filename_crop_DAPI)
    }

    if(annotation=="on"){
      print("from the file:")
      print(file_dna)
      display(img_orig)
      print("I cropped this cell:")
      display(new_img)
      print("using this mask")
      display(tmp_img)
      print("whose cell number is")
      print(cell_count)
    }
    #### strand related stuff here
  },
  error = function(e) {
    if(annotation=="on"){
      #str(e) # #prints structure of exception
      print("couldn't crop it since cell is on the edge. Neglected the following mask of a cell candidate:")
      display(tmp_img)
    }
  }
  )
}


#' get_blobs
#'
#' Makes mask of all objects bright enough
#'
#' @param img_orig Original image
#'
#' @param blob_factor, Contrast factor to multiply original image by before smoothing/smudging
#' @param bg_blob_factor, Contrast factor to multiply original image by to take background. Used prior to thresholding.
#' @param offset, Pixel value offset from bg_blob_factor. Used in thresholding to make blob mask.
#' @param final_blob_amp, Contrast factor to multiply smoothed/smudged image. Used in thresholding to make blob mask.
#' @param brush_size_blob, Brush size for smudging the dna channel to make blobs
#' @param sigma_blob, Sigma in Gaussian brush for smudging the dna channel to make blobs
#' @return Mask with cell candidates

get_blobs <- function(img_orig, blob_factor, bg_blob_factor, offset,final_blob_amp, brush_size_blob,sigma_blob){
  thresh <- blob_factor*img_orig
  # subfunction: big blur to blobs
  img_tmp_dna <- img_orig
  img_tmp <- thresh
  w <- makeBrush(size = brush_size_blob, shape = 'gaussian', sigma = sigma_blob)
  img_flo <- filter2(img_tmp, w)
  ## default amplification
  bg <- mean(bg_blob_factor*img_tmp)
  blob_th <- final_blob_amp*img_flo > bg + offset
  return(blob_th)
}


#' keep_cells
#'
#' Deletes objects in mask which are too small, large, oblong i.e. unlikely to be a cell
#'
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
    pixel_area <- x$s.area[counter]
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

