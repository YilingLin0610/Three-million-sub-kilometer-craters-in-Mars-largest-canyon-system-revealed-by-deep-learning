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
from buffer_split import buffer_split
from diagonal_lines_generation import diagonal_lines_generation
from curve_lines_check import extract_lines_with_conditions
from concurrent.futures import ProcessPoolExecutor,as_completed
from polor_projection import reproject_raster_to_polar,reproject_shp_to_polar
from tqdm import tqdm
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
    SFS_result_path=SFS_path+"result_data"
    meta_file_path=SFS_path+"test_meta"
    splits=mosaic_ID.split("_")
    latitude=int(splits[-1][1:])
    merge_shapefile_path=SFS_path+splits[0]+"_CTX_"+splits[2]+"_"+splits[3]+"_"+splits[4]+"_SeamMap_merged.shp"

    SAM_shapefile_orinigal = SFS_path+mosaic_ID+".shp"
    # # # Read the original shapefiles and split
    split_shapefilename = SAM_shapefile_orinigal[0:-4] + "_buffer_split.shp"
    if(os.path.exists(split_shapefilename)):
        pass
    else:
        split_shapefilename=buffer_split(SAM_shapefile_orinigal, SFS_path)
    # Filter out the features smaller than 7850 & features with roundness smaller than 0.83
    gdf = gpd.read_file(split_shapefilename)
    #print("lens",len(gdf))
    gdf['Area'] = gdf['geometry'].area
    tif_file = [SFS_path + x for x in os.listdir(SFS_path) if x.endswith("resample.tif")][0]

    if latitude>60 or latitude<-60: #Polor region
        #Reproject raster
        if(os.path.exists(tif_file[0:-4]+"_pro.tif")):
            pass
        else:
            reproject_raster_to_polar(
                in_tif=tif_file,
                out_tif=tif_file[0:-4] + "_pro.tif",
                pole=latitude,
                lon0=0.0,
                lat_ts=None,  # 走 Variant A（k=1）
                k=1.0
            )


        reproject_shp_to_polar(SFS_path+mosaic_ID+".shp",
                                     SFS_path+mosaic_ID+"_pro.shp", pole=latitude, lon0=0.0)
        filtered_gdf = gdf[(gdf['Area'] >= 7850) & (gdf['Area'] < 3140000)]
        filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_larger_polygon.shp", driver='ESRI Shapefile')
         #Transform to polor projection
        reproject_shp_to_polar(SAM_shapefile_orinigal[0:-4] + "_larger_polygon.shp",
                                     SAM_shapefile_orinigal[0:-4] + "_larger_polygon_pro.shp", pole=latitude, lon0=0.0)
        #Filter by circularity
        gdf_larger_polygon_pro=gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_larger_polygon_pro.shp")
        gdf_larger_polygon_pro['circularity'] = gdf_larger_polygon_pro['geometry'].apply(calculate_circularity)
        filter_gdf_larger_polygon_pro=gdf_larger_polygon_pro[gdf_larger_polygon_pro['circularity']>= roundness_thrshold]
        filter_gdf_larger_polygon_pro.to_file(SAM_shapefile_orinigal[0:-4] + "_large_pro.shp")
        filter_gdf_larger_polygon = filtered_gdf[
            (gdf_larger_polygon_pro['circularity'] >= roundness_thrshold).to_numpy()
        ]
        #Reproject to original projection for SFS processing
        filter_gdf_larger_polygon.to_file(SAM_shapefile_orinigal[0:-4] + "_large.shp")

        filtered_gdf=gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_large.shp")
    else:
        gdf['circularity'] = gdf['geometry'].apply(calculate_circularity)
        filtered_gdf = gdf[
            (gdf['Area'] >= 7850) & (gdf['circularity'] >= roundness_thrshold) & (gdf['Area'] < 3140000)]
        filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_large.shp", driver='ESRI Shapefile')





    #print("done")
    #~~filtered_gdf = gdf[(gdf['Area'] >= 7850) & (gdf['circularity'] >= roundness_thrshold)&(gdf['Area'] < 3140000)]
    #print("done2")
    #filtered_gdf=gdf

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
    filtered_gdf["fid_copy"] = filtered_gdf.index
    print(filtered_gdf.keys())

    filtered_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp", driver='ESRI Shapefile')
    print(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp")

    #filtered_gdf= gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_filter_large.shp")



    SFS_results=os.listdir(SFS_result_path)
    print("...................."+str(len(SFS_results)))




    gdf = gpd.read_file(merge_shapefile_path)
    crs=gdf.crs
    merged_gdf = gpd.GeoDataFrame()
    gdfs=[]
    gdfs_0p4=[]
    gdfs_0p6 = []
    gdfs_0p7 = []
    if(len(filtered_gdf)>0):
        for index, row in tqdm(gdf.iterrows()):
            print(row.PRODUCT_ID)
            name = row.PRODUCT_ID
            print(name)
            # Read the corresponding meta file
            print(meta_file_path + "/" + name + ".mat")
            if os.path.exists(meta_file_path + "/" + name + ".mat"):
                # print("?")
                meta_path = meta_file_path + "/" + name + ".mat"
                meta_data = scipy.io.loadmat(meta_path)
                meta_data = meta_data[name][0][0]
                azimuth_angle = meta_data[1][0][0]
                meta = meta_data[2][0][0]
                # print(meta)
                width = meta[3][0][0]
                height = meta[4][0][0]

                merge_results = np.zeros((height, width))
                # print(height,width)
                # print(name)
                SFS_name_path = []
                for i in range(len(SFS_results)):
                    # print("number",i)
                    if (name in SFS_results[i]):
                        # print("in",i)
                        SFS_name_path.append(SFS_result_path + '/' + SFS_results[i])

                # SFS_name_path = [SFS_result_path + '/' + s for s in SFS_results if name in s]
                # print(len(SFS_name_path))
                # print(len(SFS_name_path))
                for SFS in SFS_name_path:
                    #print(SFS)
                    data = scipy.io.loadmat(SFS)
                    # print(data)
                    results = data['results']
                    row_start = data['batchinfo'][0][0][1][0][0]
                    row_end = data['batchinfo'][0][0][2][0][0]
                    col_start = data['batchinfo'][0][0][3][0][0]
                    col_end = data['batchinfo'][0][0][4][0][0]
                    # print(row_start,row_end,col_start,col_end)
                    merge_results[row_start:row_end, col_start:col_end] = results
                # print("done")
                # plt.imshow(merge_results)
                # plt.show()

                row_gdf = gpd.GeoDataFrame([row], geometry='geometry')

                row_gdf.set_crs(filtered_gdf.crs, inplace=True)
                # print(row_gdf.crs)
                # print(filtered_gdf.crs)
                # #row_gdf.to_crs(filtered_gdf.crs)
                clipped_gdf = gpd.overlay(row_gdf, filtered_gdf, how='intersection')
                # print(clipped_gdf)

                # print(len(clipped_gdf))
                # clipped_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large_"+str(index), driver='ESRI Shapefile')
                # Create diagonal lines
                # print(type(SAM_shapefile_orinigal))
                name = SAM_shapefile_orinigal.split("/")

                # print(name)
                # print(SAM_shapefile_orinigal[0:-4] + "_filter_large_"+str(index)+"/"+name[-1][0:-4]+".shp")
                # #clipped_gdf = gpd.read_file(SAM_shapefile_orinigal[0:-4] + "_filter_large_" + str(index) + "/" + name[-1][
                #                                                                                                  0:-4] + "_filter_large_" + str(
                #     index) + ".shp")
                clipped_gdf['area'] = clipped_gdf.geometry.area
                clipped_gdf = clipped_gdf[clipped_gdf['area'] >= 1000]
                # print(".1231231231231212313213213")
                #print(clipped_gdf["fid_copy"])
                lines_gdf = diagonal_lines_generation(clipped_gdf, azimuth_angle)
                #print(lines_gdf["fid_copy"])
                # lines_gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large_" + str(index) + "_line.shp",
                #                   driver='ESRI Shapefile')
                # reproject_shp_to_polar(SAM_shapefile_orinigal[0:-4] + "_filter_large_" + str(index) + "_line.shp",
                #                        SAM_shapefile_orinigal[0:-4] + "_filter_large_" + str(index) + "_line_pro.shp",
                #                        pole=latitude, lon0=0.0)
                # print(type(meta))
                # print(meta[-1])
                # 提取前 6 个参数
                meta = meta[-1]
                a, b, c = meta[0, :3]  # 第一行
                d, e, f = meta[1, :3]  # 第二行

                # 创建 Affine 对象
                affine_obj = Affine(a, b, c, d, e, f)
                extracted_list = extract_lines_with_conditions(lines_gdf, merge_results, affine_obj, curve_threshold)
                # extracted_list_0p6 = extract_lines_with_conditions(lines_gdf, merge_results, affine_obj, 0.6)
                # extracted_list_0p7 = extract_lines_with_conditions(lines_gdf, merge_results, affine_obj, 0.7)
                gdf = clipped_gdf.loc[extracted_list]
                # gdf_0p6 = clipped_gdf.loc[extracted_list_0p6]
                # gdf_0p7 = clipped_gdf.loc[extracted_list_0p7]
                # 获取矩阵的形状
                height, width = merge_results.shape
                print(SAM_shapefile_orinigal[0:-4])
                s = os.path.dirname(SAM_shapefile_orinigal)
                try:
                    os.mkdir(s + "/tifs/")
                except:
                    pass
                output_tif = s + "/tifs/" + row.PRODUCT_ID + "_SFS.tif"

                # 创建 GeoTIFF 文件
                with rasterio.open(
                        output_tif,
                        'w',  # 写入模式
                        driver='GTiff',  # 驱动类型
                        height=height,  # 行数
                        width=width,  # 列数
                        count=1,  # 波段数
                        dtype=merge_results.dtype,  # 数据类型
                        crs=crs,  # 坐标系
                        transform=affine_obj  # 仿射变换矩阵
                ) as dst:
                    # 写入数据
                    dst.write(merge_results, 1)  # 将 matrix 写入第一个波段
                reproject_raster_to_polar(
                    in_tif=output_tif,
                    out_tif=output_tif[0:-4] + "_pro.tif",
                    pole=latitude,
                    lon0=0.0,
                    lat_ts=None,  # 走 Variant A（k=1）
                    k=1.0
                )
                # print(output_tif)
                # merged_gdf = gpd.pd.concat([gdf1, gdf2, gdf3], ignore_index=True)
                # gdf = gdf.loc[extracted_list]

                # gdf.to_file(SAM_shapefile_orinigal[0:-4] + "_filter_large_" + str(index) + "/" + name[-1][
                #                                                                                  0:-4] + "_filter_large_" + str(
                #     index) + "_result_update.shp", driver='ESRI Shapefile')
                # plt.imshow(merge_results)
                # plt.show()
                print(len(gdf))
                gdfs.append(gdf)
                # gdfs_0p6.append(gdf_0p6)
                # gdfs_0p7.append(gdf_0p7)






    merged_gdf = pd.concat(gdfs, ignore_index=True)
    merged_gdf.to_file(SFS_path+"merged_result_update"+"_"+str(curve_threshold)+"_roundness_"+str(roundness_thrshold)+".shp")
    if(abs(latitude)>60):
        reproject_shp_to_polar(SFS_path+"merged_result_update"+"_"+str(curve_threshold)+"_roundness_"+str(roundness_thrshold)+".shp",
                           SFS_path+"merged_result_update"+"_"+str(curve_threshold)+"_roundness_"+str(roundness_thrshold)+"_pro.shp", latitude, lon0=0.0)








