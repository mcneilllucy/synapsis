#' auto_crop_fast
#'
#' crop an image around each viable cell candidate.
#'
#' This function takes all images in a directory, and crops around individual
#' cells according to the antibody that stains synaptonemal complexes e.g.
#' SYCP3. First, it increases the brightness and smudges the image with a Gaussian
#' brush, and creates a mask using thresholding (get_blobs).
#' Then it deletes cell candidates in the mask deemed too large, too small,
#' or too long (keep_cells). Using the computeFeatures functions from
#' EBImage to locate centre and radius, the cropping area is determined and
#' the original image cropped. These images are saved in either a user specified
#' directory, or a crops folder at the location of the image files.
#'
#' @importFrom stats median sd
#' @importFrom EBImage bwlabel channel colorLabels computeFeatures
#' computeFeatures.basic computeFeatures.moment computeFeatures.shape
#' display filter2 makeBrush readImage rgbImage rmObjects rotate writeImage
#' distmap
#' watershed resize
#' @importFrom graphics text
#' @importFrom utils str
#' @export auto_crop_fast
#' @param img_path, path containing image data to analyse
#' @param path_out, user specified output path. Defaults to img_path
#' @param max_cell_area, Maximum pixel area of a cell candidate
#' @param min_cell_area, Minimum pixel area of a cell candidate
#' @param mean_pix, Mean pixel intensity of cell crop (in SYCP3 channel)
#' for normalisation
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param blob_factor, Contrast factor to multiply original image by before
#' smoothing/smudging
#' @param bg_blob_factor, Contrast factor to multiply original image by to
#' take background. Used prior to thresholding.
#' @param offset, Pixel value offset from bg_blob_factor. Used in thresholding
#' to make blob mask.
#' @param final_blob_amp, Contrast factor to multiply smoothed/smudged image.
#'  Used in thresholding to make blob mask.
#' @param test_amount, Optional number of first N images you want to run
#' function on. For troubleshooting/testing/variable calibration purposes.
#' @param brush_size_blob, Brush size for smudging the synaptonemal complex channel to make blobs
#' @param cell_aspect_ratio Maximum aspect ratio of blob to be defined as a cell
#' @param sigma_blob, Sigma in Gaussian brush for smudging the synaptonemal complex channel
#' to make blobs
#' @param channel1_string String appended to the files showing the channel
#' illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel
#' illuminating synaptonemal complexes. Defaults to SYCP3
#' @param channel3_string Optional. String appended to the files showing the
#' channel illuminating cell structures. Defaults to DAPI, if
#' third channel == "on".
#' @param third_channel Optional, defaults to "off". Set to "on" if you would
#' also like crops of the third channel.
#' @param strand_amp multiplication of strand channel for get_blobs function.
#' @param file_ext file extension of your images e.g. tif jpeg or png.
#' @param resize_l length for resized image
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in strand channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to
#' 0.05.
#' @param crowded_cells TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have many cells in a frame that almost touch
#' @param cropping_factor size of cropping window square, as factor of
#' characteristic blob radius. Defaults to 1. May need to increase if using
#' watershed.
#' @examples demo_path = paste0(system.file("extdata",package = "synapsis"))
#' auto_crop_fast(demo_path, annotation = "on", max_cell_area = 30000,
#' min_cell_area = 7000, file_ext = "tif",crowded_cells = TRUE)
#' @author Lucy McNeill
#' @return cropped synaptonemal complex and foci channels around single cells, regardless of stage

auto_crop_fast <- function(img_path,  max_cell_area = 20000, min_cell_area = 7000, mean_pix = 0.08, annotation = "off", blob_factor = 15, bg_blob_factor = 10,  offset = 0.2, final_blob_amp = 10, test_amount = 0,brush_size_blob = 51, sigma_blob = 15, channel3_string = "DAPI", channel2_string = "SYCP3", channel1_string = "MLH3", file_ext = "jpeg", third_channel = "off",cell_aspect_ratio = 2, strand_amp = 2, path_out = img_path, resize_l = 720, crowded_cells = "FALSE", watershed_radius = 50, watershed_tol = 0.2, cropping_factor =1.3)
{
  file_list <- list.files(img_path)
  dir.create(paste0(path_out,"/crops"))
  dir.create(paste0(path_out,"/crops-RGB"))
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  antibody3_store <- 0
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
      file_sc <- img_file
      print(file)
      image <- readImage(file_sc)
      img_orig <- channel(image, "grey")
      img_orig_highres <- img_orig
      img_orig <- resize(img_orig, w = resize_l, h = resize_l)
      antibody1_store <- 1
      if(annotation == "on"){
        cat("\n displaying enhanced image of cell channel")
        plot(rgbImage(blob_factor*img_orig, 0*img_orig, 0*img_orig))
      }
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
      blob_th <- get_blobs(strand_amp*img_orig,blob_factor, bg_blob_factor, offset,final_blob_amp,brush_size_blob, sigma_blob,  watershed_tol, watershed_radius, crowded_cells, annotation)
      if(crowded_cells != "TRUE"){
        blob_label <- bwlabel(blob_th)
        blob_label <- channel(blob_label, "gray")
        candidate <- bwlabel(blob_label)
      }
      else{
        candidate <- blob_th
      }

      ## remove things that aren't cells
      retained <- keep_cells(candidate, max_cell_area, min_cell_area,cell_aspect_ratio, crowded_cells, annotation)
      ### crop over each cell
      ## Loop over all objects in this final "removed" image
      if(crowded_cells == "TRUE"){
        ## do watershed on the mask again
        y <- distmap(retained)
        retained <- watershed(y)
        ##
        x_final <- computeFeatures.shape(retained)
        moment_info <- computeFeatures.moment(retained,as.matrix(img_orig))
        x_final <- data.frame(x_final)
        moment_info <- as.data.frame(moment_info)
        OOI_final <- nrow(x_final)
        counter_final <- 0
        r_max <- x_final$s.radius.max
        cx <- moment_info$m.cx
        cy <- moment_info$m.cy
      }
      else{
        x_final <- computeFeatures.shape(retained)
        moment_info <- computeFeatures.moment(retained,as.matrix(img_orig))
        x_final <- data.frame(x_final)
        moment_info <- as.data.frame(moment_info)
        OOI_final <- nrow(x_final)
        counter_final <- 0
        r_max <- x_final$s.radius.max
        cx <- moment_info$m.cx
        cy <- moment_info$m.cy
      }
      # looping through each object to crop
      while(counter_final<OOI_final){
        counter_final <- counter_final+1
        ### row of interest is the counter_final'th row of x_final
        cell_count <- cell_count +1
        if(third_channel=="on"){
          crop_single_object_fast(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI,file_sc,file_foci,file_DAPI,cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel,path_out, img_orig_highres, resize_l,crowded_cells, cropping_factor)
        }
        else{
          crop_single_object_fast(retained,OOI_final,counter_final,img_orig,img_orig_foci,img_orig_foci,file_sc,file_foci,file_foci,cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel,path_out, img_orig_highres, resize_l,crowded_cells, cropping_factor)
        }
      }
      antibody1_store <- 0
      antibody2_store <- 0
      antibody3_store <- 0
    }
  }
crop_count <- nrow(as.data.frame(list.files(paste0(img_path,"/crops-RGB/"))))
print(crop_count)
cat("out of",image_count,"images, we got",crop_count,"viable cells \n", sep = " ")
}


#' crop_single_object_fast
#'
#' Creates mask for every individual cell candidate in mask
#'
#' @param retained Mask of cell candidates which meet size criteria.
#' After smoothing/smudging and thresholding.
#' @param OOI_final, Objects of interest count. Total number of cell
#' candidates in retained.
#' @param counter_final, Counter for single cell we are focussing on. Remove
#' all other cells where counter_single not equal to counter_final.
#' @param img_orig, description
#' @param img_orig_foci, description
#' @param img_orig_DAPI, description
#' @param file_sc, filename of synaptonemal complex channel image
#' @param file_foci, filename of foci channel image
#' @param file_DAPI, filename of DAPI channel image
#' @param cell_count, counter for successful crops around cells
#' @param mean_pix, Mean pixel intensity of cell crop (in SYCP3 channel)
#' for normalisation
#' @param r_max maximum radius of blob for cropping
#' @param cx centre of blob x
#' @param cy centre of blob y
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' @param file_base, filename base common to all three channels
#' i.e. without -MLH3.jpeg etc.
#' @param img_path, path containing image data to analyse
#' @param path_out, user specified output path. Defaults to img_path
#' @param channel1_string String appended to the files showing the channel
#' illuminating foci. Defaults to MLH3
#' @param channel2_string String appended to the files showing the channel
#' illuminating synaptonemal complexes. Defaults to SYCP3
#' @param channel3_string Optional. String appended to the files showing the
#' channel illuminating cell structures. Defaults to DAPI,
#' if third channel == "on".
#' @param third_channel Optional, defaults to "off". Set to "on" if you would
#' also like crops of the third channel.
#' @param img_orig_highres the original strand image with original resolution
#' @param file_ext file extension of your images e.g. tif jpeg or png.
#' @param resize_l length of square to resize original image to.
#' @param crowded_cells TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have many cells in a frame that almost touch
#' @param cropping_factor size of cropping window square, as factor of
#' characteristic blob radius. Defaults to 1. May need to increase if using
#' watershed.
#' @return Crops around all candidates in both channels
#'
crop_single_object_fast <- function(retained, OOI_final,counter_final,img_orig,img_orig_foci,img_orig_DAPI="blank",file_sc,file_foci,file_DAPI = "blank",cell_count, mean_pix, annotation, file_base, img_path, r_max, cx, cy,channel3_string,channel2_string,channel1_string,file_ext,third_channel,path_out, img_orig_highres, resize_l, crowded_cells,cropping_factor){
  tmp_img <- retained
  ###added
  y <- distmap(tmp_img)
  blob_th <- watershed(y)
  candidate <- blob_th
  ### end added
  ## have a single object
  ### delete all other objects
  counter_single <- 0
  # looping over all other objects to crop
  if(crowded_cells == "TRUE"){
    while(counter_single < OOI_final){
      counter_single <- counter_single + 1
      # iteratively remove all other objects
      if(counter_single != counter_final){
        retained <- as.numeric(retained)*rmObjects(candidate, counter_single, reenumerate = TRUE)
        #tmp_img <- as.numeric(tmp_img)*rmObjects(bwlabel(retained), counter_single, reenumerate = TRUE)
      }
    }
    tmp_img <- retained
  }
  else{
    while(counter_single < OOI_final){
      counter_single <- counter_single + 1
      # iteratively remove all other objects
      if(counter_single != counter_final){
        tmp_img <- as.numeric(tmp_img)*rmObjects(bwlabel(retained), counter_single, reenumerate = TRUE)
      }
    }
  }

  #### resize your images here?
  dim_orig <- dim(img_orig_highres)
  new_l <- as.integer(dim_orig[1])
  noise_gone_highres <- bwlabel(resize(tmp_img, h =  new_l, w =  new_l))*as.matrix(img_orig_highres)
  noise_gone_foci <- bwlabel(resize(tmp_img, h =  new_l, w =  new_l))*as.matrix(img_orig_foci)
  if(third_channel == "on"){
    noise_gone_DAPI <- bwlabel(resize(tmp_img, h =  new_l, w =  new_l))*as.matrix(img_orig_DAPI)
  }
  crop_r_highres <- round(cropping_factor*floor(r_max[counter_final]*round( new_l/resize_l)))
  cx_highres <- cx[counter_final]*( new_l/resize_l)
  cy_highres <- cy[counter_final]*( new_l/resize_l)
  top_left_x_highres <- floor(cx_highres-crop_r_highres)
  top_left_y_highres <- floor(cy_highres-crop_r_highres)
  bottom_left_x_highres <- floor(cx_highres-crop_r_highres)
  bottom_left_y_highres <- floor(cy_highres+crop_r_highres)
  bottom_right_x_highres <- floor(cx_highres+crop_r_highres)
  bottom_right_y_highres <- floor(cy_highres+crop_r_highres)
  top_right_x_highres <- floor(cx_highres-crop_r_highres)
  top_right_y_highres <- floor(cy_highres+crop_r_highres)
  ## crop image
  ix <- bottom_left_x_highres:bottom_right_x_highres
  iy <- top_left_y_highres:bottom_left_y_highres
  ####### end high res stuff
  #########

### cropping finished
  ## cropping part
  tryCatch({
    img_path_out <- path_out
    new_img <- noise_gone_highres[ix, iy]
    ## want all images to have the same mean (mean_pix)
    orig_mean <- mean(new_img)
    mean_factor <- mean_pix/orig_mean
    new_img <- new_img*mean_factor
    file_stub <- paste0('-',channel2_string,'.',file_ext)
    file_sc <- gsub(file_stub,'', file_base)
    filename_crop <- paste0(img_path_out,"/crops/", file_sc,"-crop-",cell_count,file_stub)
    writeImage(new_img, filename_crop)
    new_img_foci <- noise_gone_foci[ix, iy]
    file_stub <- paste0('-',channel1_string,'.',file_ext)
    file_foci <- gsub(file_stub,'', file_foci)
    filename_crop_foci <- paste0(img_path_out,"/crops/", file_sc,"-crop-",cell_count,file_stub)
    writeImage(new_img_foci, filename_crop_foci)
    ### add RGB channel
    ch1 <-channel(new_img,"grey")
    ch2 <- channel(new_img_foci,"grey")
    RGB_img <- rgbImage(ch1,ch2,0*ch1)
    filename_crop_RGB <- paste0(img_path_out,"/crops-RGB/", file_sc,"-crop-",cell_count,file_stub)
    writeImage(RGB_img, filename_crop_RGB)
    if(third_channel == "on"){
      new_img_DAPI <- noise_gone_DAPI[ix, iy]
      file_stub <- paste0('-',channel3_string,'.',file_ext)
      file_DAPI <- gsub(file_stub,'', file_DAPI)
      filename_crop_DAPI <- paste0(img_path_out,"/crops/", file_sc,"-crop-",cell_count,file_stub)
      writeImage(new_img_DAPI, filename_crop_DAPI)
    }

    if(annotation=="on"){
      print("from the file:")
      print(file_sc)
      plot(img_orig)
      print("I cropped this cell:")
      plot(new_img)
      print("using this mask")
      plot(tmp_img)
      print("whose cell number is")
      print(cell_count)
    }
    #### strand related stuff here
  },
  error = function(e) {
    if(annotation=="on"){
      #str(e) # #prints structure of exception
      print("couldn't crop it since cell is on the edge. Neglected the following mask of a cell candidate:")
      plot(tmp_img)
      plot(noise_gone_highres)
    }
  }
  )
}


#' get_blobs
#'
#' Makes mask of all objects bright enough to be classified as a cell
#' cadidate
#'
#'
#'
#' @param img_orig Original image
#'
#' @param blob_factor, Contrast factor to multiply original image by before
#' smoothing/smudging
#' @param bg_blob_factor, Contrast factor to multiply original image by to take
#' background. Used prior to thresholding.
#' @param offset, Pixel value offset from bg_blob_factor. Used in thresholding
#' to make blob mask.
#' @param final_blob_amp, Contrast factor to multiply smoothed/smudged image.
#' Used in thresholding to make blob mask.
#' @param brush_size_blob, Brush size for smudging the synaptonemal complex channel to make blobs
#' @param sigma_blob, Sigma in Gaussian brush for smudging the synaptonemal complex channel to
#' make blobs
#' @param watershed_radius Radius (ext variable) in watershed method used
#' in strand channel. Defaults to 1 (small)
#' @param watershed_tol Intensity tolerance for watershed method. Defaults to 0.05.
#' @param crowded_cells TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' have many cells in a frame that almost touch
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' have many cells in a frame that almost touch
#' @return Mask with cell candidates

get_blobs <- function(img_orig, blob_factor, bg_blob_factor, offset,final_blob_amp, brush_size_blob,sigma_blob, watershed_tol, watershed_radius, crowded_cells,annotation){
  thresh <- blob_factor*img_orig
  # subfunction: big blur to blobs
  img_tmp_sc <- img_orig
  img_tmp <- thresh
  w <- makeBrush(size = brush_size_blob, shape = 'gaussian', sigma = sigma_blob)
  img_flo <- filter2(img_tmp, w)
  ## default amplification
  bg <- mean(bg_blob_factor*img_tmp)
  blob_th <- final_blob_amp*img_flo > bg + offset

  if(crowded_cells != "TRUE"){
    blob_th <- final_blob_amp*img_flo > bg + offset

  }
  else{
    y <- distmap(blob_th)
    blob_th <- watershed(y)
  }
  if(annotation == "on"){
    cat("\n here is the mask")
    if(crowded_cells == "TRUE"){
      tryCatch({
        if(nrow(data.frame(blob_th)) > 0){
          plot(colorLabels(blob_th))
        }
      },
      error = function(e) {
      }
      )
    }
    else{
      plot(blob_th)
    }

  }

  return(blob_th)
}


#' keep_cells
#'
#' Deletes objects in mask which are too small, large, oblong
#' i.e. unlikely to be a cell
#'
#' @param candidate Mask of individual cell candidates
#' @param max_cell_area, Maximum pixel area of a cell candidate
#' @param min_cell_area, Minimum pixel area of a cell candidate
#' @param cell_aspect_ratio Maximum aspect ratio of blob to be defined as a cell
#' @param crowded_cells TRUE or FALSE, defaults to FALSE. Set to TRUE if you
#' @param annotation, Choice to output pipeline choices (recommended to knit)
#' have many cells in a frame that almost touch

#' @return Mask of cell candidates which meet size criteria
keep_cells <- function(candidate, max_cell_area, min_cell_area, cell_aspect_ratio, crowded_cells, annotation){

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
  if(annotation == "on"){
    cat("\n displaying the retained cells in mask (correct size/ shape)")
    if(crowded_cells == "TRUE"){
      tryCatch({
        if(nrow(data.frame(retained)) > 0){
          plot(retained)
        }
      },
      error = function(e) {
      }
      )

    }
    else{
      plot(retained)
    }
  }
  retained <- bwlabel(retained)
  return(retained)
}

