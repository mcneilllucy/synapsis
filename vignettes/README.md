Welcome to the synapsis vignette homepage.

Here is a short guide to navigate the vignettes with a flowchart questionnaire.

![questionnaire](resources/figures/questionnaire.png)

Generally, your analysis will involve (1) data preparation, (2) auto cropping, meiosis stage identification (optional) and (3) foci counting and/or measuring distances between foci.

For (1), please go through the yellow part of the flowchart with reference to `setting_up_synapsis.Rmd`.

For (2), please go through the magenta part of the flowchart with reference to `synapsis_tutorial_v1.Rmd`.

If your cell is type A, you can call

```r
auto_crop_fast(path, annotation = "on", max_cell_area = 30000, min_cell_area = 7000)
```

Otherwise if your cell is type B, try something like

```r
auto_crop_fast(path, annotation = "on", max_cell_area = 300000, min_cell_area = 7000, channel3_string = "c1", channel2_string = "c3", channel1_string = "c2", file_ext = "tif",cell_aspect_ratio = 5)
```

For (3) please go through the blue part of the flowchart and keep looking at `synapsis_tutorial_v1.Rmd`.

if you are looking at mouse pachytene (what synapsis is designed for), the following crop or distance functions can be called with these parameters:

```r
SYCP3_stats <- get_pachytene(path,ecc_thresh = 0.8, area_thresh = 0.04)
foci_counts <- count_foci(path,offset_factor = 5, brush_size = 1, brush_sigma = 1, annotation = "on",WT_out = "Fancm+/+",KO_out = "Fancm-/-", stage = "pachytene")
```

and

```r
df_dist <- measure_distances_general(path,offset_factor = 5, brush_size = 1, brush_sigma = 1, annotation = "on", WT_out = "Fancm+/+",KO_out = "Fancm-/-", stage = "pachytene")
```

otherwise, for everything else, you could call something like:

```r
foci_counts <- count_foci(path,offset_factor = 5, brush_size = 1, brush_sigma = 1, annotation = "on",WT_str = "tif",WT_out ="NA",file_ext = "tif", channel2_string = "c3", channel1_string = "c2", offset_px = 0.4, crowded_foci = TRUE)
```
and

```r
df_dist <- measure_distances_general(path,offset_factor = 5, brush_size = 1, brush_sigma = 1, annotation = "on",WT_str = "tif",WT_out ="NA",file_ext = "tif", channel2_string = "c3", channel1_string = "c2", offset_px = 0.4, crowded_foci = TRUE)
```
