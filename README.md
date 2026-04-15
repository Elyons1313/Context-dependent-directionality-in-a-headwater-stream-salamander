2025_Combined_Position_Data.csv - full capture history dataset, without including habitat information
  #this sheet was used when developing the model, but adding the habitat data changed the subsetting and I elected to use this one to filter the data instead of rewriting the entire thing

2025_With_Habitat - same thing, just with Mesohabitat data - only used for Mesohabitat data 

2025_Capture_History.csv - earlier version of the dataset, not used in this paper

Abundance_JS_Bear.R - The code for the Open Population Jolly Seber Model used in the analyses
  _Excluding_Telemetry - same model, but for paradise brook
  _Zigzag.R - same model, but for Zigzag Brook

HBFish_Data_Combined.csv - the data for all the brook trout captured from 2022-2025. 

Publication_Code.RMD - R markdown document containing all the publication code (outside of the abundance model) and the figures used in the paper 
  # the abundance model was excluded because it takes a long time to run (~a few hours) and wanted to keep the publication code easily runnable

Relative_Abundances_2025.csv - output from the abundance model, used to allow the publication code to run quickly. Combined the outputs from all three streams into one csv
  # running one of the individual abundance models should give the same (ish) output as this csv for each stream - confirm if desired

Stream_Data_Winsor - all the wetted width data measurements, combined into one neat csv without summary statistics. (Winsor is my advisor's name, sent it to him while I was still in the field and didn't rename it)
