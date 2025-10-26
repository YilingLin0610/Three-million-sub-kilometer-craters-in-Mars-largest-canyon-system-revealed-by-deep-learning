import scipy.io
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from shapely.geometry import mapping
from osgeo import gdal
import os
import math
import numpy as np
from affine import Affine
import pandas as pd
import rtree
import matplotlib.pyplot as plt
import glob
from buffer_split_60 import buffer_split
from diagonal_lines_generation import diagonal_lines_generation
from curve_lines import extract_lines_with_conditions

"""
0. Process the SAM output TIFF files
  0.1 Transform the TIFF files to shapefile format
  0.2 Split the original shapefile files (Consider whether this step smooth the features)
  0.3 Delete the features whose areas smaller than 7850m² and whose roundness smaller than 0.83
1. Read the SFS results and merge to different scenes
2. Clip the shapefile to the extent of different scenes
3. filter out through SFS results
"""
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

def clip_with_rtree(filtered_gdf, clip_feature):

    # 创建 RTree 索引
    idx = rtree.index.Index()
    for i, geom in enumerate(filtered_gdf.geometry):

        idx.insert(i, geom.bounds)
    #print("done2")

    # 使用 RTree 索引查询可能与裁剪范围相交的几何对象
    possible_matches_index = list(idx.intersection(clip_feature))
    #print("done3")
    possible_matches = filtered_gdf.iloc[possible_matches_index]
    #print(len(possible_matches))

    # # 对筛选出的几何对象进行裁剪
    # # 修复自相交
    # possible_matches =possible_matches.buffer(0)
    # clipped_gdf = possible_matches.clip(clip_feature)
    #
    # # 返回裁剪后的结果和几何对象数量
    # return clipped_gdf, len(clipped_gdf)








def post_process(mosaic_ID, SFS_path):
    SFS_result_path=SFS_path+"result_data_60"

    splits=mosaic_ID.split("_")
    merge_shapefile_path=SFS_path+splits[0]+"_CTX_"+splits[2]+"_"+splits[3]+"_"+splits[4]+"_SeamMap_merged_buffer.shp"
    #print(merge_shapefile_path)
    #SAM_shapefile_orinigal=SFS_path+splits[0]+"_CTX_"+splits[2]+"_"+splits[3]+"_"+splits[4]+"_Mosaic_resample.shp"
    SAM_shapefile_orinigal = SFS_path+"results_60/merge_prediction_shp/"+mosaic_ID+".shp"
    print(SAM_shapefile_orinigal)
    # Read the original shapefiles and split
    #split_shapefilename = buffer_split(SAM_shapefile_orinigal, SFS_path)
    split_shapefilename=SAM_shapefile_orinigal[0:-4] + "_buffer_split.shp"
    #gdf.to_file(split_shapefilename)
    gdf=gpd.read_file(SAM_shapefile_orinigal)
    gdf=gdf[(gdf['Area'] >= 196250)]
    gdf=gdf.buffer(-200)
    gdf=gdf.buffer(200)
    gdf.to_file(split_shapefilename)
    gdf = gpd.read_file(split_shapefilename)
    #print("lens",len(gdf))
    gdf['Area'] = gdf['geometry'].area

    gdf['circularity'] = gdf['geometry'].apply(calculate_circularity)
    #print("done")
    filtered_gdf = gdf[(gdf['Area'] >= 502400) & (gdf['circularity'] >= 0.85)&(gdf['Area'] < 314000000)]
    #print("done2")
    filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4]+"_large.shp", driver='ESRI Shapefile')
    #print(SAM_shapefile_orinigal[0:-4]+"_large")
    #print("done3")
    #filtered_gdf_line = gdf_line[filter_condition]

    # 用于存储需要保留的要素索引
    to_keep = []
    deleted=[]

    # 创建空间索引
    spatial_index = filtered_gdf.sindex

    # 遍历所有要素
    for i, geom in enumerate(filtered_gdf.geometry):
        #print(i)
        contained = False
        if geom is None or geom is None:
            continue  # 或者可以选择标记为需要删除：to_keep.add(i)
            # 使用空间索引查询可能相交的要素
        possible_matches_index = list(spatial_index.intersection(geom.bounds))
        #print("possible_matches_index",possible_matches_index)
        #possible_matches = filtered_gdf.iloc[possible_matches_index]

            # 遍历可能相交的要素
        for j in possible_matches_index:
            other_geom=filtered_gdf.geometry.iloc[j]
            #print("list possible",j)
            if i != j:
                if(j in deleted):
                        pass
                else:
                    try:
                        intersection_area = geom.intersection(other_geom).area
                    except:
                        geom = geom.buffer(0)
                        other_geom = other_geom.buffer(0)
                        intersection_area = geom.intersection(other_geom).area
                    #print(intersection_area / geom.area)
                    # 如果重叠面积超过当前要素面积的 95%，标记其他要素为需要删除
                    if intersection_area >= 0.85 * geom.area:
                        contained = True  # 保留当前要素
                        deleted.append(i)
                        break  # 跳出内层循环

        to_keep.append(not contained)
        # 如果当前要素已经被标记为需要删除，跳过



    #print(to_keep)
    filtered_gdf = filtered_gdf.iloc[list(to_keep)]

    filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp", driver='ESRI Shapefile')
    filtered_gdf= gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp")
    #
    #
    #
    # SFS_results=os.listdir(SFS_result_path)
    # print(SFS_results)
    # print("...................."+str(len(SFS_results)))
    #crs = filtered_gdf.crs
    # for i in range(len(SFS_results)):
    #
    #     output_tif=r"F:\Mars_SFS\CTX_images\MurrayLab_GlobalCTXMosaic_V01_E-088_N-04\tifs_60\\"+name+".tif"
    #     with rasterio.open(
    #             output_tif,
    #             'w',  # 写入模式
    #             driver='GTiff',  # 驱动类型
    #             height=height,  # 行数
    #             width=width,  # 列数
    #             count=1,  # 波段数
    #             dtype=matrix.dtype,  # 数据类型
    #             crs=crs,  # 坐标系
    #             transform=affine_obj  # 仿射变换矩阵
    #     ) as dst:
    #         # 写入数据
    #         dst.write(matrix, 1)  # 将 matrix 写入第一个波段






    gdf = gpd.read_file(merge_shapefile_path)

    crs=gdf.crs
    gdfs=[]

    for index,row in gdf.iterrows():
        #print(row)
        name=row['PRODUCT_ID']

        data = scipy.io.loadmat(SFS_result_path + "/" + name+".mat")
        data = data[name][0][0]
        azimuth_angle = data[1][0]
        #print(azimuth_angle)
        # print(azimuth_angle)
        matrix = data[4]
        #print(np.shape(matrix))
        # print(data[2])
        meta = data[2][0][0]
        # print(meta)
        width = meta[3][0][0]
        height = meta[4][0][0]
        meta = meta[-1]
        a, b, c = meta[0, :3]  # 第一行
        d, e, f = meta[1, :3]  # 第二行

        # 创建 Affine 对象
        affine_obj = Affine(a, b, c, d, e, f)

        row_gdf = gpd.GeoDataFrame([row], geometry='geometry')
        row_gdf.set_crs(filtered_gdf.crs, inplace=True)
        # print(row_gdf.crs)
        # print(filtered_gdf.crs)
        # #row_gdf.to_crs(filtered_gdf.crs)
        clipped_gdf = gpd.overlay(row_gdf, filtered_gdf, how='intersection')
        #clipped_gdf.to_file(SFS_path+"merged_result_update_60_1.shp")
        # print(len(clipped_gdf))
        # clipped_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large_"+str(index), driver='ESRI Shapefile')
        # Create diagonal lines
        # print(type(SAM_shapefile_orinigal))
        name = SAM_shapefile_orinigal.split("/")

        lines_gdf = diagonal_lines_generation(clipped_gdf, azimuth_angle)

        extracted_list = extract_lines_with_conditions(lines_gdf, matrix, affine_obj, 0.5)
        # extracted_list_0p6 = extract_lines_with_conditions(lines_gdf, merge_results, affine_obj, 0.6)
        # extracted_list_0p7 = extract_lines_with_conditions(lines_gdf, merge_results, affine_obj, 0.7)
        gdf = clipped_gdf.loc[extracted_list]
        gdfs.append(gdf)






    merged_gdf = pd.concat(gdfs, ignore_index=True)
    print(len(merged_gdf))


    # 用于存储需要保留的要素索引
    to_keep = []
    deleted=[]

    # 创建空间索引
    spatial_index = merged_gdf.sindex

    # 遍历所有要素
    for i, geom in enumerate(merged_gdf.geometry):
        #print(i)
        contained = False
        if geom is None or geom is None:
            continue  # 或者可以选择标记为需要删除：to_keep.add(i)
            # 使用空间索引查询可能相交的要素
        possible_matches_index = list(spatial_index.intersection(geom.bounds))
        #print("possible_matches_index",possible_matches_index)
        #possible_matches = filtered_gdf.iloc[possible_matches_index]

            # 遍历可能相交的要素
        for j in possible_matches_index:
            other_geom=merged_gdf.geometry.iloc[j]
            #print("list possible",j)
            if i != j:
                if(j in deleted):
                        pass
                else:
                    try:
                        intersection_area = geom.intersection(other_geom).area
                    except:
                        geom = geom.buffer(0)
                        other_geom = other_geom.buffer(0)
                        intersection_area = geom.intersection(other_geom).area
                    #print(intersection_area / geom.area)
                    # 如果重叠面积超过当前要素面积的 95%，标记其他要素为需要删除
                    if intersection_area >= 0.85 * geom.area:
                        contained = True  # 保留当前要素
                        deleted.append(i)
                        break  # 跳出内层循环

        to_keep.append(not contained)
    merged_gdf = merged_gdf.iloc[list(to_keep)]
        # 如果当前要素已经被标记为需要删除，跳过


    merged_gdf.to_file(SFS_path+"merged_result_update_60_0.85.shp")
    print("saved")

