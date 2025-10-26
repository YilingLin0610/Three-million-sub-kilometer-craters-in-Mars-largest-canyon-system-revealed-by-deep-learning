import os
import geopandas as gpd
from osgeo import gdal,ogr
import numpy as np

def buffer_split(shapefilename, SFS_path):
    # 打印输入文件名
    print(os.path.basename(shapefilename))

    # 读取 Shapefile
    gdf = gpd.read_file(shapefilename)
    #gdf = gdf[(gdf['Area'] <= 1500000)]
    print(len(gdf))
    #gdf = gdf.head(1000)

    # 初始化几何列表和标志列表
    geometries = gdf.geometry.tolist()  # 将几何对象转换为列表
    geometries_flag = np.zeros(len(geometries), dtype=int)  # 标志列表

    # 处理几何对象
    while np.any(geometries_flag == 0):
        # 找到第一个未处理的几何对象
        index = np.where(geometries_flag == 0)[0][0]
        #print("Processing geometry at index:", index)


        geometry = geometries[index]
        distance = -200  # 初始缓冲距离
        buffered_geometry = geometry.buffer(distance, resolution=2)

        if buffered_geometry.is_empty:
            geometries_flag[index] = 1  # 标记为已处理
            geometry_2 = geometry.buffer(-200, resolution=2)
            geometries[index] = geometry_2.buffer(200, resolution=2)
        else:
            while not buffered_geometry.is_empty:
                if buffered_geometry.geom_type == 'MultiPolygon':  # 如果缓冲后几何变为 MultiPolygon
                    # 移除原始几何对象
                    geometries.pop(index)
                    geometries_flag = np.delete(geometries_flag, index)

                    # 将拆分后的多边形添加到几何列表
                    for polygon in buffered_geometry.geoms:
                        single_polygon = polygon.buffer(-distance, resolution=2)
                        geometries.append(single_polygon)
                        geometries_flag = np.append(geometries_flag, 0)
                    break
                else:
                    # 减少缓冲距离
                    distance -= 200
                    buffered_geometry = geometry.buffer(distance, resolution=2)

                    if buffered_geometry.is_empty:
                        geometry_2 = geometry.buffer(-200, resolution=2)
                        geometries[index] = geometry_2.buffer(200, resolution=2)

                        geometries_flag[index] = 1  # 标记为已处理

    # 创建新的 GeoDataFrame
    result_gdf = gpd.GeoDataFrame(geometry=geometries, crs=gdf.crs)

    # 保存结果
    output_path = shapefilename.replace(".shp", "_buffer_split.shp")
    result_gdf.to_file(output_path)

    return output_path


