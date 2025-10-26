from merge_600 import post_process
import os
print("yes")
if __name__ == "__main__":
    total_path = r"F:\Mars_SFS\CTX_images_Valles_600\\"
    SFS_path = total_path + 'results_60\merge_prediction_shp\\'
    SAM_shapefile_orinigal=total_path + 'results_60\merge_prediction_shp\CTX_images_Valles_600.shp'
    post_process(SAM_shapefile_orinigal, SFS_path)


    #mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N-04'

    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-036_N-20'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N00'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)


    # MurrayLab_CTX_V01_E104_N - 04
    # _SeamMap_merged.shp
