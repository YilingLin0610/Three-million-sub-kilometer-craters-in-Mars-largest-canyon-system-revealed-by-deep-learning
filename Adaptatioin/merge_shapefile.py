
import os
import geopandas as gpd
import os
import pandas as pd


if __name__ == "__main__":

    total_path=r"/data/A/Yiling/Mars_SAM_crater/CTX_images_Valles//"
    Mosaic_names = [ f for f in os.listdir(total_path) if os.path.isdir(os.path.join(total_path, f))]



    # 初始化列表保存所有 GeoDataFrame
    gdf_list = []

    num=0
    for mosaic_ID in Mosaic_names:
        if(mosaic_ID=='MurrayLab_GlobalCTXMosaic_V01_E104_N-04'):
            print("skip")
        else:
            print(mosaic_ID)
            shp = total_path + mosaic_ID + "/" + "merged_result_update_60_0.85.shp"
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
    total_path=total_path+"merged_result_update_60_0.85_all.shp"
    merged_gdf.to_file(total_path)


    print(f"✅ 合并完成，共 {num} 个要素")





    #mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N-04'

    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-036_N-20'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N00'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)


    # MurrayLab_CTX_V01_E104_N - 04
    # _SeamMap_merged.shp
