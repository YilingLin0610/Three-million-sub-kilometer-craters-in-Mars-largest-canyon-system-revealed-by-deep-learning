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
from buffer_split_600 import buffer_split
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








def post_process(SAM_shapefile_orinigal,SFS_path):

    # Read the original shapefiles and split
    split_shapefilename=SAM_shapefile_orinigal[0:-4] + "_buffer_split.shp"
    #gdf=gpd.read_file(SAM_shapefile_orinigal)
    #split_shapefilename = buffer_split(SAM_shapefile_orinigal, SFS_path)
    gdf = gpd.read_file(split_shapefilename)
    gdf['Area'] = gdf['geometry'].area
    gdf['circularity'] = gdf['geometry'].apply(calculate_circularity)
    filtered_gdf = gdf[(gdf['Area'] >= 28260000) & (gdf['circularity'] >= 0.85)&(gdf['Area'] < 1256000000)]
    filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4]+"_large.shp", driver='ESRI Shapefile')


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
    lines_gdf = diagonal_lines_generation(filtered_gdf, 45)
    lines_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "buffer_result_diagonal_lines.shp")

    tiff_file = r"F:\Mars_SFS\CTX_images_Valles_600\MOLA_DEM_pro.tif"
    with rasterio.open(tiff_file) as src:
        matrix = src.read(1)  # 读取第一波段为 numpy array（从 1 开始计）
        print(matrix)
        affine_obj = src.transform  # 获取仿射变换对象 Affine
        crs = src.crs  # 获取坐标参考系统（可选）
    extracted_list = extract_lines_with_conditions(lines_gdf, matrix, affine_obj, 300)
    merged_gdf = filtered_gdf.loc[extracted_list]
    merged_gdf.to_file(SFS_path + "merged_result_update_600_f300_85.shp")



