#' get_pachytene
#'
#' Identifies crops in pachytene
#'
#' @import EBImage
#' @import stats
#' @import graphics
#' @import utils
#' @export get_pachytene
#' @param file_list The file list
#' @param img_path The path
#' @return Pairs of foci and SC channel crops for pachytene


get_pachytene <- function(img_path, species_num = 20, offset = 0.2,ecc_thresh = 0.85, area_thresh = 0.06)
{
  # input :

  # output : a bunch of output jpegs? Or save them all?



  #BiocManager::install("EBImage")
  library(EBImage)
  cell_count <- 0
  image_count <-0
  antibody1_store <- 0
  antibody2_store <- 0
  pachytene_count <- 0
  setwd(img_path)
  max_obj <- species_num + round(0.1*species_num)

  df_cols <- c("filename","cell_no","genotype", "px_mask","px_total", "px_fraction", "mean_ecc","mean_ratio","skew","sd_bright_px","stage_classification")
  df_cells <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
  colnames(df_cells) <- df_cols

  img_path_new <- paste0(img_path,"/crops/")
  print(img_path_new)

  setwd(img_path_new)
  dir.create("pachytene")
  file_list <- list.files(img_path_new)
  print(file_list)

  ## for each image that is *-dna.jpeg,
  for (file in file_list){
    setwd(img_path_new)
    if(grepl("*SYCP3.jpeg", file)){
      file_dna = file
      image_count <- image_count +1
      image <- readImage(file_dna)
      img_orig <- channel(image, "grey")
      antibody1_store <- 1
    }
    if(grepl("*MLH3.jpeg", file)){
      file_foci = file
      #print(file_foci)
      image <- readImage(file_foci)
      img_orig_foci <- channel(image, "gray")
      # call functions: get
      antibody2_store <- 1
    }
    if(antibody1_store + antibody2_store ==2){
      antibody1_store <- 0
      antibody2_store <- 0
      print("I have a pair")
      new_img<-img_orig
      #### now see which have the right amount of strands
      disc = makeBrush(21, "disc")
      disc = disc / sum(disc)
      localBackground = filter2(new_img, disc)

      thresh_crop = (new_img - localBackground > offset)
      strands <- bwlabel(thresh_crop)
      color_img_strands<- colorLabels(strands, normalize = TRUE)
      num_strands <- computeFeatures.shape(strands)
      num_strands <- data.frame(num_strands)
      print(nrow(num_strands))
      #### segment the strands
      if (nrow(num_strands)<max_obj && nrow(num_strands)>5){
        cell_count <- cell_count + 1
        ### identified a good image. count foci
        #display(new_img)
        #display(strands)

        ###### adding data frame info


        ### data frame stuff

        if(grepl( "++", file, fixed = TRUE) == TRUE){
          genotype <- "Fancm+/+"
        }

        if(grepl( "--", file, fixed = TRUE) == TRUE){
          genotype <- "Fancm-/-"
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
        #moment_info <- computeFeatures.haralick(bwlabel(tmp_img),as.matrix(noise_gone))
        moment_info <- as.data.frame(moment_info)
        ## might actually want to find the real centre first..




        ### data frame stuff ends
        dim_img <- nrow(new_img)
        #dim_img <- dim_img[1,1]
        num_px <- dim_img*dim_img
        px_mask <- sum(SC_areas)
        px_total <- num_px
        px_fraction <- px_mask/px_total
        mean_ecc <- mean(moment_info$m.eccentricity)
        stage_classification <- "pachytene"

        sd_bright_px <- sd(image_mat)





        ######

        if(mean_ecc > ecc_thresh){
          if(px_fraction > area_thresh){

            df_cells <- rbind(df_cells,t(c(file,cell_count,genotype,px_mask, px_total,px_fraction, mean_ecc,mean_ratio,skew,sd_bright_px,"stage_classification")))

            pachytene_count <- pachytene_count + 1

            file_dna <- tools::file_path_sans_ext(file_dna)
            filename_crop = paste0("./pachytene/", file_dna,".jpeg")
            writeImage(img_orig, filename_crop)

            file_foci <- tools::file_path_sans_ext(file_foci)
            filename_crop_foci = paste0("./pachytene/", file_foci, ".jpeg")
            writeImage(img_orig_foci, filename_crop_foci)

            print("two filenames are")
            print(filename_crop)
            print(filename_crop_foci)
            print("end filename")

          }
        }



      }
      ###
    }
  }
print("number of cells kept")
print(pachytene_count)
colnames(df_cells) <- df_cols
return(df_cells)
}

