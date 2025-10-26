
import os
import geopandas as gpd
import os
import pandas as pd


if __name__ == "__main__":

    total_path=r"F:\Mars_SFS\CTX_images//"
    files=[r"F:\Mars_SFS\CTX_images\merged_result_update_60_all_rectify_merge_6_filter.shp",
           r"F:\Mars_SFS\CTX_images_Valles_600\results_60\merge_prediction_shp\merged_result_update_600_f300_85_rectified.shp",
           r"F:\Mars_SFS\CTX_images\MurrayLab_GlobalCTXMosaic_V01_E-040_N-20\results_60\merge_prediction_shp\MurrayLab_GlobalCTXMosaic_V01_E-040_N-20_filter_large.shp",
           r"F:\Mars_SFS\CTX_images\MurrayLab_GlobalCTXMosaic_V01_E-056_N-16\results_60\merge_prediction_shp\MurrayLab_GlobalCTXMosaic_V01_E-056_N-16_filter_large.shp"

    ]



    # 初始化列表保存所有 GeoDataFrame
    gdf_list = []

    num=0
    for mosaic_ID in files:
        if(mosaic_ID=='MurrayLab_GlobalCTXMosaic_V01_E104_N-04'):
            print("skip")
        else:
            print(mosaic_ID)
            shp = mosaic_ID
            # shp = total_path + mosaic_ID + "/results_60/merge_prediction_shp/" + mosaic_ID+"_large.shp"
            # if(os.path.exists ( total_path + mosaic_ID + "/"+"merged_result_update_60_new.shp")):
            gdf = gpd.read_file(shp)
            len_gdf = len(gdf)
            num = num + len_gdf
            print(num)


        gdf_list.append(gdf)

    # 合并所有 GDF
    merged_gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)

    # 保存到新文件
    total_path=total_path+"results_LT_1km_original.shp"
    merged_gdf.to_file(total_path)


    print(f"✅ 合并完成，共 {num} 个要素")


