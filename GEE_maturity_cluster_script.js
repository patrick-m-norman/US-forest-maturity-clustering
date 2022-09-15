//The following script was used to produce the US forest maturity clustering using Google Earth Engine 
//In total, five additional datasets were required as inputs:
//GLOBAL
//Tree cover - ee.Image("UMD/hansen/global_forest_change_2020_v1_8")
//Tree height - ee.ImageCollection("users/potapovpeter/GEDI_V27")
//Above ground biomass sourced from https://data.globalforestwatch.org/datasets/gfw::aboveground-live-woody-biomass-density/about

//US
//Level 3 ecoregions (polygonized using GDAL) sourced from https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states
//FIA forest type groups sourced from https://data.fs.usda.gov/geodata/rastergateway/forest_type/

//The only additional data included was a raster template for exporting and a list of all the forest type groups per ecoregions, which is
//publicly available as the asset ee.FeatureCollection("users/patrickmnorman/Forest_and_ecoregion_for_GEE_table")

//The conus forest type group raster can be downloaded from https://data.fs.usda.gov/geodata/rastergateway/forest_type/
//----------------------------------------------------------------------------------------------------------------------------------------
//Loading in each of the datasets
var Tree_cover = ee.Image("UMD/hansen/global_forest_change_2020_v1_8"),
    Tree_height = ee.ImageCollection("users/potapovpeter/GEDI_V27"),
    Above_ground_biomass = ee.Image("users/carlycampbell/gbf_biomass"),
    L3_ecoregions = ee.FeatureCollection("users/patrickmnorman/us_eco_l3_projected"),
    forest_groups = ee.Image("users/patrickmnorman/conus_forest_type_group_raster"),
    forest_group_by_ecoregion_table = ee.FeatureCollection("users/patrickmnorman/Forest_and_ecoregion_for_GEE_table"),
    export_area = 
    ee.Geometry.Polygon(
        [[[-120.34269150604206, 49.411586310089774],
          [-125.17667588104206, 49.124835932771155],
          [-126.05558213104206, 42.42415213066956],
          [-120.91398056854206, 31.87835785573225],
          [-112.74015244354206, 31.166615920449363],
          [-105.66495713104206, 29.38299619625266],
          [-100.78702744354206, 26.70720166805698],
          [-97.00773056854206, 24.887291326544982],
          [-87.33976181854206, 24.2478237627109],
          [-79.73722275604206, 24.687808623264196],
          [-79.51749619354206, 29.574276646005014],
          [-80.52823838104206, 30.563072375440527],
          [-76.61710556854206, 33.688565850869765],
          [-74.33194931854206, 36.103137824059104],
          [-73.54093369354206, 39.13079122682691],
          [-69.19034775604206, 40.94742542972901],
          [-67.52042588104206, 43.03746462571629],
          [-66.15812119354206, 44.55983485421651],
          [-66.55362900604206, 46.710382098450566],
          [-67.47648056854206, 47.9611334169127],
          [-70.46476181854206, 47.93169773498136],
          [-74.37589463104206, 46.19569439857067],
          [-77.62784775604206, 44.65369451468965],
          [-81.71476181854206, 42.68312813059799],
          [-81.53898056854206, 44.965464671473796],
          [-82.98917588104206, 46.589716680606614],
          [-85.88956650604206, 48.37147370192815],
          [-91.03116806854206, 49.00966892455189],
          [-95.90909775604206, 49.86692413351389],
          [-102.94034775604206, 49.69667136125509]]]);



//Loading in each of the indices used for clustering

var forest_cover_projections = Tree_cover.select(['treecover2000']);
var biomass = Above_ground_biomass.select(['b1']);
var height = Tree_height.toBands().select(['GEDI_NAM_v27_b1']);

//Getting the projection and resolution details about the cover dataset to then merge the three datasets
var CoverProjection = forest_cover_projections.projection();
var CoverRes = forest_cover_projections.projection().nominalScale().getInfo()

//Aligning the three datasets to tree cover
var forest_groupRescaled = forest_groups
      // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mode(),
      maxPixels: 65535
    })
    // Request the data at the scale and projection of the height image
    .reproject({
      crs: CoverProjection,
      scale: CoverRes
    });


var BiomassRescaled = biomass
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })
    .reproject({
      crs: CoverProjection,
      scale: CoverRes
    });
    
var HeightRescaled = height
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 65535
    })
    .reproject({
      crs: CoverProjection,
      scale: CoverRes
    });
    


//The function below loops through the forest by ecoregion table supplied .
var Maturity_cluster = forest_group_by_ecoregion_table.aggregate_array('ogc_fid').sort().slice(1,873,1)/*list_subset if needed*/.map(function(n) {
  //Accessing the row then ecoregion and forest type group names
  var table_row = ee.Feature(forest_group_by_ecoregion_table.filterMetadata('ogc_fid', 'equals', n).first())
  var forest_type = ee.Number(table_row.get('dn'))
  var ecoregion = table_row.get('ecoregion_')
  var ecoregion_poly = L3_ecoregions.filterMetadata('US_L3NAME', 'equals', ecoregion);
  //Once the forest type group has been identified, it is then clipped to the ecoregion
  var forest_group_clip = forest_groupRescaled.clip(ecoregion_poly)
    
    //Masking each of the input metrics to the forest type group by ecoregion
    var forest_cover_clip = forest_cover_projections.updateMask(forest_group_clip.eq(forest_type)).clip(ecoregion_poly);
    var BiomassRescaled_clip = BiomassRescaled.updateMask(forest_group_clip.eq(forest_type)).clip(ecoregion_poly);
    var HeightRescaled_clip = HeightRescaled.updateMask(forest_group_clip.eq(forest_type)).clip(ecoregion_poly);
    
    //Generating the 0th, 25th, 50th, 75th and 100th percentiles for tree cover, tree heigh and above groud biomass layers
    var CoverPercentiles = forest_cover_clip.reduceRegion({
      reducer: ee.Reducer.percentile([0,25,50,75,100]), 
      geometry: ecoregion_poly,
      scale: CoverRes,
      maxPixels: 1e13
    });

    var BiomassPercentiles = BiomassRescaled_clip.reduceRegion({
      reducer: ee.Reducer.percentile([0,25,50,75,100]),
      geometry: ecoregion_poly,
      scale: CoverRes,
      maxPixels: 1e13
      });
    var heightPercentiles = HeightRescaled_clip.reduceRegion({
      reducer: ee.Reducer.percentile([0,25,50,75,100]),
      geometry: ecoregion_poly,
      scale: CoverRes,
      maxPixels: 1e13
      });
      
    //Providing a score for each of the quartiles. 1 = below 25%, 2 = between 25% and 50%, 3 = between 50% and 75% and 4 = above 75%
    var Forest_cover_classq = forest_cover_clip
                            .where(forest_cover_clip.gte(ee.Number(CoverPercentiles.get('treecover2000_p0'))).and(forest_cover_clip.lt(ee.Number(CoverPercentiles.get('treecover2000_p25')))), 1)
                            .where(forest_cover_clip.gte(ee.Number(CoverPercentiles.get('treecover2000_p25'))).and(forest_cover_clip.lt(ee.Number(CoverPercentiles.get('treecover2000_p50')))), 2)
                            .where(forest_cover_clip.gte(ee.Number(CoverPercentiles.get('treecover2000_p50'))).and(forest_cover_clip.lt(ee.Number(CoverPercentiles.get('treecover2000_p75')))), 3)
                            .where(forest_cover_clip.gte(ee.Number(CoverPercentiles.get('treecover2000_p75'))), 4)
                            
    var Biomas_classq = BiomassRescaled_clip
                            .where(BiomassRescaled_clip.gte(ee.Number(BiomassPercentiles.get('b1_p0'))).and(BiomassRescaled_clip.lt(ee.Number(BiomassPercentiles.get('b1_p25')))), 1)
                            .where(BiomassRescaled_clip.gte(ee.Number(BiomassPercentiles.get('b1_p25'))).and(BiomassRescaled_clip.lt(ee.Number(BiomassPercentiles.get('b1_p50')))), 2)
                            .where(BiomassRescaled_clip.gte(ee.Number(BiomassPercentiles.get('b1_p50'))).and(BiomassRescaled_clip.lt(ee.Number(BiomassPercentiles.get('b1_p75')))), 3)
                            .where(BiomassRescaled_clip.gte(ee.Number(BiomassPercentiles.get('b1_p75'))), 4)
                            
    var Height_classq = HeightRescaled_clip
                            .where(HeightRescaled_clip.gte(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p0'))).and(HeightRescaled_clip.lt(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p25')))), 1)
                            .where(HeightRescaled_clip.gte(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p25'))).and(HeightRescaled_clip.lt(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p50')))), 2)
                            .where(HeightRescaled_clip.gte(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p50'))).and(HeightRescaled_clip.lt(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p75')))), 3)
                            .where(HeightRescaled_clip.gte(ee.Number(heightPercentiles.get('GEDI_NAM_v27_b1_p75'))), 4)                       
  
    var output_layer = Forest_cover_classq.add(Biomas_classq).add(Height_classq);
    return output_layer
});

//Converting the image to an image collection
var collection = ee.ImageCollection(Maturity_cluster).toBands()
var combined = collection.reduce(ee.Reducer.mode()).toByte()

//Mapping colours to explore output
var vis = {min: 3, max: 12, palette: 'black,white'};
Map.addLayer(combined, vis, '')


//Exporting the image for further processing/
Export.image.toDrive({
  image: combined,
  description: 'Maturity_cluster_image',
  region: export_area,
  scale:30,
  maxPixels: 1e13
});//
