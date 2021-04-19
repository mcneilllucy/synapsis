crop <- function(file_list, img_path)
{
  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  pair <- 0
  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    setwd(img_path)
    if(grepl("*dna.jpeg$", file)){
      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(image, "grey")
      pair <- 0
    }
    if(grepl("*foci.jpeg$", file)){
      file_foci = file
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      pair <- 1
    }
    if(pair ==1){

      w = makeBrush(size = 101, shape = 'gaussian', sigma = 51)
      img_flo = filter2(img_orig, w)
      disc = makeBrush(51, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(img_orig, w)
      offset = 0.2
      nucBadThresh = (img_orig - localBackground > offset)
      img_tmp_dna <- img_orig
      img_tmp <- nucBadThresh
      w = makeBrush(size = 51, shape = 'gaussian', sigma = 15)
      img_flo = filter2(img_tmp, w)
      bg <- mean(10*img_tmp)
      offset = 0.2
      blob_th = 10*img_flo > bg + offset
      blob_label = bwlabel(blob_th)
      blob_label <- channel(blob_label, "gray")
      candidate <- bwlabel(blob_label)

      # delete everything that's too small
      colorimg<- colorLabels(candidate, normalize = TRUE)
      x <- computeFeatures.shape(candidate)
      x <- data.frame(x)
      OOI <- width(x)
      counter <- 0
      removed <- candidate

      while(counter<OOI){
        counter <- counter+1
        pixel_area = x$s.area[counter]
        semi_maj <- x$s.radius.max[counter]
        semi_min <- x$s.radius.min[counter]
        # if statement checking if it's the wrong area
        if(pixel_area> 20000 | pixel_area < 5000){
          removed <- as.numeric(removed)*rmObjects(candidate, counter, reenumerate = TRUE)
        }
        ## if statement checking that it's not too long i.e. not at edge.
        if(semi_maj/semi_min > 2 & is.na(semi_maj/semi_min)==FALSE){
          removed <- as.numeric(removed)*rmObjects(candidate, counter, reenumerate = TRUE)
        }

      }
      removed <- bwlabel(removed)
      ### crop foci channel here
      ## Loop over all objects in this final "removed" image
      x_final <- computeFeatures.shape(removed)
      x_final <- data.frame(x_final)
      OOI_final <- width(x_final)
      counter_final <- 0
      # looping through each object to crop
      while(counter_final<OOI_final){
        counter_final <- counter_final+1
        tmp_img <- removed
        ## have a single object
        ### delete all other objects
        counter_single <- 0
        # looping over all other objects to crop
        while(counter_single < OOI_final){
          counter_single <- counter_single + 1
          # iteratively remove all other objects
          if(counter_single != counter_final){
            tmp_img <- as.numeric(tmp_img)*rmObjects(bwlabel(removed), counter_single, reenumerate = TRUE)

          }

        }
        noise_gone <- bwlabel(tmp_img)*as.matrix(img_orig)
        noise_gone_foci <- bwlabel(tmp_img)*as.matrix(img_orig_foci)
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
          tryCatch({
            new_img <- noise_gone[ix, iy]
            ## want all images to have the same mean to 0.1
            orig_mean <- mean(new_img)
            mean_factor <- 0.08/orig_mean
            new_img <- new_img*mean_factor
            cell_count <- cell_count + 1
            display(new_img)
            print("cell count for above crop is")
            print(cell_count)
            if(cell_count != 0){
              print("keeping it")
            }
            else{
              print("not keeping that one")
            }


            #### strand related stuff here
          },
          error = function(e) {
            #what should be done in case of exception?
            str(e) # #prints structure of exception

          }
          )
        }
      }
    }
    pair <- 0
  }

print("out of")
print(image_count)
print("we got")
print(cell_count)
print("viable cells")
}
