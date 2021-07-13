test_that("Creates a crops folder", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"")
  auto_crop_fast(demo_path)
  folders <- list.dirs(demo_path)
  for (folder in folders){
    if(grepl("crops", folder)){
      folder_exists <- TRUE
    }
  }
  expect_true(folder_exists)
})

