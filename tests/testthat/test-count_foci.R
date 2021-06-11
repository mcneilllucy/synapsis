test_that("outputs data frame which has correct number of columns", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"/")
  test_df <- count_foci(demo_path)
  expect_s3_class(test_df, "data.frame")
})

test_that("outputs data frame which has correct number of rows", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"/")
  demo_path_2 <-paste0(system.file("extdata",package = "synapsis"),"/crops/pachytene/")
  files <- list.files(demo_path_2)
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
  expect_equal(nrow(test_df), file_count/2)
})
