import geopandas as gpd
from shapely.geometry import box
from tqdm import tqdm
import numpy as np
# 读取两个 shapefile
gdf2 = gpd.read_file(r"F:\Mars_SFS\CTX_images\results_LT_1km.shp")
gdf1 = gpd.read_file(r"F:\Mars_SFS\Robbins_2012\robbins_12_valles_new_clip_pro.shp")
area_threshold = np.pi * 5000 * 5000  # ≈ 78,539,816 平方米

# 筛选出面积小于该阈值的要素
gdf1 = gdf1[gdf2.geometry.area <= area_threshold]

# 确保使用相同的坐标系
gdf1 = gdf1.to_crs(gdf2.crs)

# 创建 shapefile2 的空间索引
sindex = gdf2.sindex

# 创建一个布尔掩码，初始化为全 True（表示全部保留）
mask = [True] * len(gdf2)

# 遍历 shapefile1 的每个要素
for geom in tqdm(gdf1.geometry):
    # 获取可能与之相交的 gdf2 要素索引（通过 bounding box）
    possible_idx = list(sindex.intersection(geom.bounds))

    # 精确判断并更新掩码
    for idx in possible_idx:
        if gdf2.geometry.iloc[idx].intersects(geom):
            mask[idx] = False  # 有交集的要素标记为 False（排除）

# 筛选无交集的 gdf2 要素
gdf2_no_intersection = gdf2[mask]

# 保存结果
gdf2_no_intersection.to_file("shapefile2_no_overlap.shp")

# 保存结果（可选）
gdf2_no_intersection.to_file(r"F:\Mars_SFS\Robbins_2012\newly_found.shp")