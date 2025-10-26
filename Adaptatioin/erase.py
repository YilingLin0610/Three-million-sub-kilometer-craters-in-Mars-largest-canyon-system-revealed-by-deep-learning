import geopandas as gpd
from shapely.ops import unary_union
import os
import pandas as pd
# # 读取大要素（被擦除的基础要素）
# large_gdf = gpd.read_file(r"E:\Mars_results_0906\mxd\result_data\Valles_extent.shp")
#
# total_path = r"F:\Mars_SFS\CTX_images\\"
# Mosaic_names = [f for f in os.listdir(total_path) if os.path.isdir(os.path.join(total_path, f))]
#
# # # 初始化列表保存所有 GeoDataFrame
# # gdf_list = []
# #
# # num = 0
# # for mosaic_ID in Mosaic_names:
# #     if (mosaic_ID == 'MurrayLab_GlobalCTXMosaic_V01_E104_N-04'):
# #         print("skip")
# #     else:
# #         print(mosaic_ID)
# #         name=[f for f in os.listdir( total_path + mosaic_ID) if f.endswith("_merged.shp")]
# #         shp = total_path + mosaic_ID + "/" + name[0]
# #         # shp = total_path + mosaic_ID + "/results_60/merge_prediction_shp/" + mosaic_ID+"_large.shp"
# #         # if(os.path.exists( total_path + mosaic_ID + "/"+"merged_result_update_60_new.shp")):
# #         gdf = gpd.read_file(shp)
# #         len_gdf = len(gdf)
# #         num = num + len_gdf
# #         print(num)
# #
# #     gdf_list.append(gdf)
# #
# # # 合并所有 GDF
# # merged_gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)
# #
# # # 保存到新文件
# # total_path=total_path+"merged_source_extent.shp"
# # merged_gdf.to_file(total_path)
# merged_gdf=gpd.read_file(r"E:\Mars_results_0906\mxd\result_data\merge_extent_CTX.shp")
#
#
# #print(f"✅ 合并完成，共 {len(num)} 个要素")
#
# # 将多个小要素合并为一个几何对象
# union_small = unary_union(merged_gdf.geometry)
#
# # 从大要素中擦除这些区域
# large_gdf['geometry'] = large_gdf.geometry.apply(lambda geom: geom.difference(union_small))
#
# # 保存擦除后的结果为新的shapefile
# large_gdf.to_file(total_path+"no_Data_extent_1.shp")
import rasterio
from rasterio.features import shapes
import geopandas as gpd
from tqdm import tqdm
import numpy as np
from shapely.geometry import shape
total_path = r"F:\Mars_SFS\CTX_images\\"
Mosaic_names = [f for f in os.listdir(total_path) if os.path.isdir(os.path.join(total_path, f))]
# 输入和输出路径
tif_path = total_path + Mosaic_names[0] + "/" +[f for f in os.listdir( total_path + Mosaic_names[0]) if f.endswith("resample.tif")][0]
output_shp = total_path+"nodata_extent.shp"
area_threshold = 20000000  # ㎡
Mosaic_names=Mosaic_names
# 读取 TIF 文件
# 打开 TIF 文件
all_geometries = []
for mosaic_ID in tqdm(Mosaic_names):
    name=[f for f in os.listdir( total_path + mosaic_ID) if f.endswith("_resample.tif")]
    tif_path = total_path + mosaic_ID + "/" + name[0]
    with rasterio.open(tif_path) as src:
        data = src.read(1)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs

        if nodata is None:
            print(f"跳过无 nodata 值的文件: {tif_path}")
            continue

        mask = data == nodata
        shape_generator = shapes(data, mask=mask, transform=transform)

        for geom, value in shape_generator:
            if value == nodata:
                all_geometries.append({
                    "geometry": shape(geom),
                    "crs": crs
                })


    # 设置 CRS（使用最后一个 tif 的）
gdf = gpd.GeoDataFrame(geometry=[g["geometry"] for g in all_geometries])
gdf.set_crs(crs, inplace=True)



    # 计算面积并筛选
gdf["area_m2"] = gdf.geometry.area
gdf = gdf[gdf["area_m2"] > area_threshold]

    # 保存为合并的 Shapefile
gdf.to_file(output_shp)
print(f"处理完成，输出保存为：{output_shp}")

