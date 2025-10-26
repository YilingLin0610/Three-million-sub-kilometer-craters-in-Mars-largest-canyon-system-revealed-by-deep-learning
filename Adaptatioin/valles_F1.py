import os
import geopandas as gpd
import geopandas as gpd
import numpy as np
import cv2
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import os
from scipy.optimize import leastsq


# 圆拟合函数：最小二乘法
def fit_circle_least_squares(x, y):
    def calc_R(xc, yc):
        return np.sqrt((x - xc) ** 2 + (y - yc) ** 2)

    def f(c):
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    center_estimate = (np.mean(x), np.mean(y))
    center, _ = leastsq(f, center_estimate)
    xc, yc = center
    Ri = calc_R(xc, yc)
    R = Ri.mean()
    return xc, yc, R




def count_null_production_id(directory):
    null_count = 0  # 记录 null 值的数量
    total_count = 0
    total_count_original=0

    # 遍历文件夹内所有 shapefile 文件
    for filename in os.listdir(directory):
        if filename.endswith("_rectified.shp"):  # 确保文件名以 "_rectified" 结尾
            filepath = os.path.join(directory, filename)
            # 读取 shapefile 文件
            gdf = gpd.read_file(filepath)
            total_count+=len(gdf)
            print(total_count)

            # 统计 'PRODUCTION_id' 字段为 null 的要素数量
            null_count += gdf['PRODUCT_ID'].isnull().sum()
            print(null_count)
        if filename.endswith(".shp") and not filename.endswith("_rectified.shp"):
            filepath = os.path.join(directory, filename)
            gdf = gpd.read_file(filepath)
            total_count_original+=len(gdf)
            print(total_count_original)

    return null_count

# 输入主目录路径
Total_path = r"F:\Mars_SFS\Samples_test_valles/results"

# 遍历子文件夹
Mosaic_names = [
    name for name in os.listdir(Total_path)
    if os.path.isdir(os.path.join(Total_path, name))
]

for Mosaic_name in tqdm(Mosaic_names):
    print(f"处理：{Mosaic_name}")
    Filepath = os.path.join(Total_path, Mosaic_name)
    shp_path = os.path.join(Filepath, "merged_result_update.shp")

    if not os.path.exists(shp_path):
        print(f"未找到 shapefile：{shp_path}")
        continue

    gdf = gpd.read_file(shp_path)

    # 面积过滤（排除太大的）
    threshold_area = 500 * 500 * 3.14
    gdf = gdf[gdf.geometry.area < threshold_area]

    results = []
    circles = []

    for idx, geom in enumerate(gdf.geometry):
        if geom.geom_type == 'Polygon':
            coords = np.array(geom.exterior.coords)
        elif geom.geom_type == 'MultiPolygon':
            coords = np.array(list(geom.geoms)[0].exterior.coords)
        else:
            continue

        if len(coords) < 5:
            continue

        try:
            x_coords = coords[:, 0]
            y_coords = coords[:, 1]
            x_center, y_center, radius = fit_circle_least_squares(x_coords, y_coords)
        except Exception as e:
            print(f"拟合失败（索引 {idx}）：{e}")
            continue

        # 拟合圆边界点
        theta = np.linspace(0, 2 * np.pi, 100)
        circle_x = x_center + radius * np.cos(theta)
        circle_y = y_center + radius * np.sin(theta)
        circle_points = np.column_stack((circle_x, circle_y))
        circle_polygon = Polygon(circle_points)

        # 记录属性
        results.append({
            "index": idx,
            "center_x": x_center,
            "center_y": y_center,
            "radius": radius
        })

        circles.append({
            "geometry": circle_polygon,
            "center_x": x_center,
            "center_y": y_center,
            "radius": radius
        })

    # 保存为 GeoDataFrame（Shapefile）
    if circles:
        gdf_circles = gpd.GeoDataFrame(circles)
        gdf_circles.set_crs(gdf.crs, inplace=True)
        gdf_circles.to_file(os.path.join(Filepath, "circle_fits.shp"))

    # 保存为 CSV
    if results:
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(Filepath, "circle_fits.csv"), index=False)
# 示例文件夹路径
directory = r"F:\Mars_SFS\Samples_test_valles/results"  # 请替换为实际文件夹路径
result = count_null_production_id(directory)
print(f"Total null 'PRODUCTION_id' elements: {result}")