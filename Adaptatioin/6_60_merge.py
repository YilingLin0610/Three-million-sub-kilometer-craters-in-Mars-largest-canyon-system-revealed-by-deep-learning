import geopandas as gpd
import pandas as pd
import os
import numpy as np
total_path = r"F:\Mars_SFS\CTX_images\\"
Mosaic_names = [f for f in os.listdir(total_path) if os.path.isdir(os.path.join(total_path, f))]
def calculate_circularity(geometry):
    """
    计算几何对象的圆度。
    :param geometry: Shapely 几何对象
    :return: 圆度值
    """
    #print(geometry)
    if geometry is None:
        return np.nan
    area = geometry.area  # 计算面积
    perimeter = geometry.length  # 计算周长
    if perimeter == 0:  # 防止除以零
        return np.nan
    return (4 * np.pi * area) / (perimeter ** 2)
# 初始化列表保存所有 GeoDataFrame
gdf_list = []

num = 0
# for mosaic_ID in Mosaic_names:
#     if (mosaic_ID == 'MurrayLab_GlobalCTXMosaic_V01_E104_N-04'):
#         print("skip")
#     else:
#         print(mosaic_ID)
#         shp = total_path + mosaic_ID + "/" + "merged_result_update.shp"
#         # shp = total_path + mosaic_ID + "/results_60/merge_prediction_shp/" + mosaic_ID+"_large.shp"
#         # if(os.path.exists( total_path + mosaic_ID + "/"+"merged_result_update_60_new.shp")):
#         gdf = gpd.read_file(shp)
#         len_gdf = len(gdf)
#         num = num + len_gdf
#         print(num)
#
#     gdf_list.append(gdf)
#
# # 合并所有 GDF
# gdf1= gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)

# 读取两个 Shapefile 文件
#gdf1 = gpd.read_file(shapefile1_path)

gdf2 = gpd.read_file(r"F:\Mars_SFS\CTX_images\merged_result_update_60_all_rectify.shp")

# 设定面积阈值
threshold_area = 500 * 500 * 3.14
#gdf1_filtered = gdf1[gdf1.geometry.area > threshold_area]
#gdf1_filtered.to_file(r"F:\Mars_SFS\CTX_images\merged_6_1km.shp")
gdf1_filtered=gpd.read_file(r"F:\Mars_SFS\CTX_images\merged_6_1km.shp")
# gdf1_filtered['circularity'] = gdf1_filtered['geometry'].apply(calculate_circularity)
#
# gdf1_filtered=gdf1_filtered['circularity'] >= 0.9
# 合并两个 GeoDataFrame
merged = gpd.GeoDataFrame(pd.concat([gdf2, gdf1_filtered], ignore_index=True), crs=gdf2.crs)

# 保存为新文件
merged.to_file(r"F:\Mars_SFS\CTX_images\merged_result_update_60_all_rectify_merge_6.shp")



# 定义函数，判断是否有超过85%的重叠
# 用于存储需要保留的要素索引
to_keep = []
deleted = []

# 创建空间索引
spatial_index = merged.sindex

# 遍历所有要素
for i, geom in enumerate(merged.geometry):
    # print(i)
    contained = False
    if geom is None or geom is None:
        continue  # 或者可以选择标记为需要删除：to_keep.add(i)
        # 使用空间索引查询可能相交的要素
    possible_matches_index = list(spatial_index.intersection(geom.bounds))
    # print("possible_matches_index",possible_matches_index)
    # possible_matches = filtered_gdf.iloc[possible_matches_index]

    # 遍历可能相交的要素
    for j in possible_matches_index:
        other_geom = merged.geometry.iloc[j]
        #print("list possible",j)
        if i != j:
            if (j in deleted):
                pass
            else:
                try:
                    intersection_area = geom.intersection(other_geom).area
                except:
                    geom = geom.buffer(0)
                    other_geom = other_geom.buffer(0)
                    intersection_area = geom.intersection(other_geom).area
                # print(intersection_area / geom.area)
                # 如果重叠面积超过当前要素面积的 95%，标记其他要素为需要删除
                if intersection_area >= 0.6 * geom.area:
                    print("yes")
                    contained = True  # 保留当前要素
                    deleted.append(i)
                    break  # 跳出内层循环

    to_keep.append(not contained)
    # 如果当前要素已经被标记为需要删除，跳过

print(len(to_keep))
print(sum(to_keep))
filtered_gdf = merged.iloc[list(to_keep)]
filtered_gdf['Area'] = filtered_gdf['geometry'].area
filtered_gdf=filtered_gdf[(filtered_gdf['Area'] >=785000)]

filtered_gdf.to_file(r"F:\Mars_SFS\CTX_images\merged_result_update_60_all_rectify_merge_6_filter.shp", driver='ESRI Shapefile')
