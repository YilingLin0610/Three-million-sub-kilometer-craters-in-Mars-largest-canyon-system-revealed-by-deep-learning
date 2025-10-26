import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import rowcol
from scipy.optimize import curve_fit
from scipy.signal import detrend
from tqdm import tqdm
import matplotlib.pyplot as plt


# 定义固定函数
def fixed_function(x, a1, b1, c1, a2, b2, c2):
    return a1 * np.exp(-((x - b1) ** 2) / (c1 )) + a2 * np.exp(-((x - b2) ** 2) / (c2 ))

def extract_pixel_value(x, y,band):
    # Transform geographic coordinates to pixel coordinates
    px = int((x - geotransform[0]) / geotransform[1])
    py = int((y - geotransform[3]) / geotransform[5])
    if(px>Xsize-1):
        px=Xsize-1
    elif(px<0):
        px=0
    if(py>Ysize-1):
        py=Ysize-1
    elif(py<0):
        py=0
    return band.ReadAsArray(px, py, 1, 1)[0, 0]


# 定义主函数
def extract_lines_with_conditions(line_gdf, raster_matrix, geotrans, thre):



    # 初始化存储提取值的列表
    line_values = []
    idxs=[]

    # 遍历每条线
    for idx, line in tqdm(line_gdf.iterrows(), total=line_gdf.shape[0]):
        #print(line.Name)
        #print(line)
        #print(idx)
        geometry = line.geometry  # 获取几何对象
        if geometry.geom_type == "LineString":
            # 提取线的点
            points = list(geometry.coords)  # 获取线的所有点坐标 [(x1, y1), (x2, y2), ...]
            #print(points)

            # 提取每个点的像素值
            values = []
            for x, y in points:
                # 将坐标转换为栅格的行列号
                row, col = rowcol(geotrans, x, y)
                if 0 <= row < raster_matrix.shape[0] and 0 <= col < raster_matrix.shape[1]:  # 检查是否在栅格范围内
                    value = raster_matrix[row, col]
                    values.append(value)
                else:
                    values.append(np.nan)  # 如果点不在栅格范围内，填充 NaN
            idxs.append(idx)

            line_values.append(values)

    # 初始化存储提取结果的列表
    extracted_list = []

    # 遍历每条线的值
    for i, values in tqdm(enumerate(line_values)):
        contains_nan = np.isnan(values).any()
        #print(values)
        if not contains_nan:
            # 去趋势
            detrended_values = detrend(values)
            detrended_values = detrended_values - np.min(detrended_values)

            # 生成 x 数据
            x_data = np.array(range(len(values)))
            x_data = x_data - np.min(x_data)
            # plt.figure()
            # plt.plot(x_data, detrended_values)

            try:
                print("yes")
                # 定义参数范围
                lower_bounds = [0, 0, 0.1, 0, np.max(x_data) / 2, 0.1]
                upper_bounds = [np.max(detrended_values), np.max(x_data) / 2, np.inf, np.max(detrended_values),
                                np.max(x_data), np.inf]

                # 拟合曲线
                params, params_covariance = curve_fit(fixed_function, x_data, detrended_values,
                                                      bounds=(lower_bounds, upper_bounds), maxfev=10000)

                # 计算拟合值
                y_fitted = fixed_function(x_data, *params)
                ss_res = np.sum((detrended_values - y_fitted) ** 2)
                ss_tot = np.sum((detrended_values - np.mean(detrended_values)) ** 2)
                r_squared = 1 - (ss_res / ss_tot)
                plt.plot(x_data, y_fitted)



                # 提取参数
                a1, b1, c1, a2, b2, c2 = params
                b_center = (b1 + b2) / 2
                y_lower = fixed_function(b_center, *params)
                height = max(a1 - y_lower, a2 - y_lower)

                b_max = max(b1, b2)
                b_min = min(b1, b2)

                # # 判断是否满足条件
                if height > thre and b_min < np.max(x_data) * 2 / 6 and b_max > np.max(x_data) * 4 / 6:
                     extracted_list.append(idxs[i])
                # if height > thre:
                #     extracted_list.append(idxs[i])
                # plt.text(0.25 * np.max(x_data), 0.8 * np.max(detrended_values), f'R² = {r_squared:.2f}')
                # plt.text(0.25 * np.max(x_data), 0.7 * np.max(detrended_values), f'height = {height:.2f}')
                # plt.text(0.25 * np.max(x_data), 0.6 * np.max(detrended_values), f'b1 = {b1:.2f}')
                # plt.text(0.25 * np.max(x_data), 0.5 * np.max(detrended_values), f'b2 = {b2:.2f}')
                # #plt.text(0.25 * np.max(x_data), 0.4 * np.max(detrended_values), f'ratio = {ratio:.2f}')
                # plt.title(str(idxs[i]))
                #plt.ylim(min(detrended_values),min(detrended_values)+1.5)
                #plt.show()
                #plt.savefig(r"E:\Mars_results_0906\CTX_images\MurrayLab_GlobalCTXMosaic_V01_E104_N-04\MurrayLab_GlobalCTXMosaic_V01_E104_N-04\MurrayLab_CTX_V01_E104_N-04_Mosaic_resample_filter_large_14\figures_idx\\"+str(idxs[i])+".jpg")

            except Exception as e:
                print(f"拟合失败: {idxs[i]}, 错误: {e}")

    # 返回提取的线索引
    return extracted_list