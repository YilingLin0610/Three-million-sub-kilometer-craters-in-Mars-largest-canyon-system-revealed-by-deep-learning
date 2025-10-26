import scipy.io
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from scipy.constants import precision
from shapely.geometry import mapping
from osgeo import gdal
import os
import math
import numpy as np
from affine import Affine
import pandas as pd
import rtree
import matplotlib.pyplot as plt
from buffer_split import buffer_split
from diagonal_lines_generation import diagonal_lines_generation
from curve_lines import extract_lines_with_conditions
from concurrent.futures import ProcessPoolExecutor,as_completed
from polor_projection import reproject_shp_to_polar,reproject_raster_to_polar
from tqdm import tqdm
import time

def affine_from_feature(feature, pixel_size_x, pixel_size_y=None):
    """
    根据 gdf 的单个 polygon 要素和已知像素大小，构造仿射变换 Affine
    feature: GeoSeries 或 shapely Polygon
    pixel_size_x: 像素宽度
    pixel_size_y: 像素高度（默认等于宽度，且为负）
    """
    geom = feature.geometry
    minx, miny, maxx, maxy = geom.bounds


    a = pixel_size_x
    e = -(pixel_size_y if pixel_size_y is not None else pixel_size_x)  # 常见情况：行号向下，取负
    b, d = 0.0, 0.0
    c, f = minx, maxy  # 左上角

    return Affine(a, b, c, d, e, f)
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








def post_process(mosaic_ID, SFS_path,curve_threshold,roundness_thrshold):
    SFS_result_path=SFS_path+"seam_post_results"
    meta_file_path=SFS_path+"test_meta"
    splits=mosaic_ID.split("_")
    latitude=int(splits[-1][1:])

    merge_shapefile_path=SFS_path+splits[0]+"_CTX_"+splits[2]+"_"+splits[3]+"_"+splits[4]+"_SeamMap_merged_post_new_buffer.shp"
    SAM_shapefile_orinigal = SFS_path+"results/merge_prediction_shp//"+mosaic_ID+".shp"
    filtered_gdf= gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp")

    SFS_results=os.listdir(SFS_result_path)
    gdf = gpd.read_file(merge_shapefile_path)

    gdfs=[]

    if(len(filtered_gdf)>0):
        for index, row in gdf.iterrows():
            parent=row.PARENT_ID
            name=parent+"_"+str(index)+".mat"
            print(name)

            if os.path.exists(SFS_result_path+"/"+name):
                #print("yes")
                #print(SFS_result_path+"/"+name)
                result_data = scipy.io.loadmat(SFS_result_path+"/"+name)
                #print(result_data)
                results=result_data['results']
                azimuth_angle=result_data['azimuth'][0][0]
                #print(np.shape(results))
                #print(azimuth_angle)
                # result_data=result_data['S'][0][0]
                # azimuth_angle=result_data[1][0][0]
                # results=result_data[0]
                # print(azimuth_angle)
                # print(np.shape(results))
                affine_obj=affine_from_feature(row, 6)

                row_gdf = gpd.GeoDataFrame([row], geometry='geometry')

                row_gdf.set_crs(filtered_gdf.crs, inplace=True)
                # print(row_gdf.crs)

                clipped_gdf = gpd.overlay(row_gdf, filtered_gdf, how='intersection')

                lines_gdf = diagonal_lines_generation(clipped_gdf, azimuth_angle)


                extracted_list = extract_lines_with_conditions(lines_gdf, results, affine_obj, curve_threshold)

                gdf = clipped_gdf.loc[extracted_list]



                print(len(gdf))
                gdfs.append(gdf)

        merged_gdf = pd.concat(gdfs, ignore_index=True)
        merged_gdf.to_file(SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
            roundness_thrshold) + "_seam_supp.shp")
    else:
        gdf_2 = gpd.GeoDataFrame(columns=["id", "name"], geometry=[], crs=gdf.crs)

        # 保存为 shapefile
        gdf_2.to_file(SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
            roundness_thrshold) + "_seam_supp.shp")







def list_dirs(path="."):
    """只列出目录，不列文件"""
    return [
        name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))
    ]
if __name__ == "__main__":
    loop = 1
    while (loop == 1):
        total_paths = [
            r"/data/C/Yiling/Mars_SFS/CTX_images_process_1/",
        ]
        for total_path in total_paths:
            print(total_path)
            for mosaic_ID in tqdm(list_dirs(total_path)):
                print(mosaic_ID)
                SFS_path = total_path + mosaic_ID + "/"
                SFS_result_path = SFS_path + "seam_post_results"
                meta_file_path_2 = SFS_path + "seam_post"
                curve_threshold = 0.5
                roundness_thrshold = 0.83
                if(os.path.exists(SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
                                roundness_thrshold) + "_final.shp")):
                    print("processed")
                else:
                    if (os.path.exists(SFS_result_path) and os.path.exists(meta_file_path_2) and os.path.exists(
                            SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
                                roundness_thrshold) + ".shp")):
                        if (len(os.listdir(SFS_result_path)) == len(os.listdir(meta_file_path_2))):
                            post_process(mosaic_ID, SFS_path, 0.5, 0.83)
                            # 读取两个shapefile

                            gdf1 = gpd.read_file(
                                SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
                                    roundness_thrshold) + "_seam_supp.shp")
                            gdf2 = gpd.read_file(
                                SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
                                    roundness_thrshold) + ".shp")

                            # 合并（行拼接）
                            filtered_gdf = gpd.GeoDataFrame(
                                pd.concat([gdf1, gdf2], ignore_index=True),
                                crs=gdf1.crs
                            )
                            #print(len(filtered_gdf))

                            # 用于存储需要保留的要素索引
                            to_keep = []
                            deleted = []

                            # 创建空间索引
                            spatial_index = filtered_gdf.sindex

                            # 遍历所有要素
                            for i, geom in enumerate(filtered_gdf.geometry):
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
                                    other_geom = filtered_gdf.geometry.iloc[j]
                                    # print("list possible",j)
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
                                            if intersection_area >= 0.85 * geom.area:
                                                contained = True  # 保留当前要素
                                                deleted.append(i)
                                                break  # 跳出内层循环

                                to_keep.append(not contained)
                            print(len(to_keep))
                            # 如果当前要素已经被标记为需要删除，跳过
                            filtered_gdf = filtered_gdf.iloc[list(to_keep)]
                            filtered_gdf.to_file(
                                SFS_path + "merged_result_update" + "_" + str(curve_threshold) + "_roundness_" + str(
                                    roundness_thrshold) + "_final.shp")






        time.sleep(600)


















