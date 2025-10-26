
import math
import geopandas as gpd
from shapely.geometry import LineString, Point

def diagonal_lines_generation(gdf,angle_deg):
    angle_rad = math.radians(angle_deg)  # 转换为弧度
    custom_index=[]
    # 遍历每个要素
    lines = []  # 用于存储所有生成的线
    fid_copy_list = []
    for idx, row in gdf.iterrows():

        #print(idx)
        geom = row.geometry  # 获取几何对象

        # 计算几何的 centroid
        centroid = geom.centroid

        # 计算几何的 bounding box
        minx, miny, maxx, maxy = geom.bounds
        width = maxx - minx  # 宽度
        height = maxy - miny  # 高度

        # 计算对角线长度的一半
        diagonal_length = math.sqrt(width ** 2 + height ** 2) * 0.5

        # 计算偏移量（相对于 Y 轴的角度）
        dx = diagonal_length * math.sin(angle_rad)  # X 偏移
        dy = diagonal_length * math.cos(angle_rad)  # Y 偏移

        # 计算线的起点和终点
        x_start = centroid.x - dx
        y_start = centroid.y - dy
        x_end = centroid.x + dx
        y_end = centroid.y + dy

        # 计算线的总长度
        line_length = math.sqrt((x_end - x_start) ** 2 + (y_end - y_start) ** 2)

        # 设置间隔（例如 6 米）
        interval = 6

        # 计算点的数量
        num_points = int(line_length // interval)  # 间隔数
        #print(num_points)

        # 创建线的点列表
        points = []
        if(num_points>0):
            for i in range(num_points + 1):  # 包括起点和终点
                # print(i)
                t = i / num_points  # 参数 t 从 0 到 1
                x = x_start + t * (x_end - x_start)
                y = y_start + t * (y_end - y_start)
                points.append((x, y))
                # print(points)
                # 将点列表转换为 LineString 对象
            line = LineString(points)
            lines.append(line)
            custom_index.append(idx)
            fid_copy_list.append(row["fid_copy"])


    # 创建 GeoDataFrame
    lines_gdf = gpd.GeoDataFrame({"fid_copy": fid_copy_list},geometry=lines, crs=gdf.crs,index=custom_index)



    return lines_gdf
