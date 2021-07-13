test_that("outputs data frame which has correct number of columns", {
  demo_path <-paste0(system.file("extdata",package = "synapsis"),"")
  test_df <- measure_distances_general(demo_path, target_foci_number = 2)
  expect_s3_class(test_df, "data.frame")
})
