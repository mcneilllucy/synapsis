test_that("Creates a pachytene folder", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"")
  auto_crop_fast(demo_path)
  SYCP3_stats <- get_pachytene(demo_path,ecc_thresh = 0.8, area_thresh = 0.04)
  folders <- list.dirs(demo_path)
  for (folder in folders){
    if(grepl("pachytene", folder)){
      folder_exists <- TRUE
    }
  }
  expect_true(folder_exists)
})

test_that("Files in pairs in pachytene folder", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"/pachytene")
  files <- list.files(demo_path)
  file_count <- 0
  for (file in files){
    if(grepl("*SYCP3.jpeg", file)){
      file_count <- file_count + 1
    }
    if(grepl("*MLH3.jpeg", file)){
      file_count <- file_count + 1
    }
  }
  test_df <- count_foci(demo_path)
  expect_equal(file_count%%2, 0)
})
