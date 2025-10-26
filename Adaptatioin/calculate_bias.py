from osgeo import gdal, ogr
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
from scipy.optimize import curve_fit
import geopandas as gpd
from tqdm import tqdm
import os
import math
def fixed_function(x, a1, b1,c1, a2,b2,c2):
    return a1 * np.exp(-((x - (b1)) ** 2)/c1)+a2 * np.exp(-((x - (b2)) ** 2)/c2)
# Define the fixed function to fit
def extract_pixel_value(x, y,band,geotransform,Xsize,Ysize):
    # Transform geographic coordinates to pixel coordinates
    # print(x,geotransform[0])
    px = int((x - geotransform[0]) / geotransform[1])
    py = int((y - geotransform[3]) / geotransform[5])
    # print( geotransform[1])
    # print(px,py)
    # print(Xsize,Ysize)

    if(px>Xsize-1):
        px=Xsize-1
    elif(px<0):
        px=0
    if(py>Ysize-1):
        py=Ysize-1
    elif(py<0):
        py=0
    # print(band.ReadAsArray(px, py, 1, 1)[0, 0])
    return band.ReadAsArray(px, py, 1, 1)[0, 0]
def calculate_circularity(geometry):
    """
    计算几何对象的圆度。
    :param geometry: Shapely 几何对象
    :return: 圆度值
    """
    area = geometry.area  # 计算面积
    perimeter = geometry.length  # 计算周长
    if perimeter == 0:  # 防止除以零
        return np.nan
    return (4 * np.pi * area) / (perimeter ** 2)
if __name__ == "__main__":


    # Input and output shapefile paths
    study_area = 8
    angle_deg = 65.341


    # Example: 45 degrees, but you can change this value
    input_shapefile = r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  ".shp"
    output_shapefile = r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_line.shp"

    # Open the input shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(input_shapefile, 0)  # Read-only
    layer = datasource.GetLayer()
    source_spatial_ref = layer.GetSpatialRef()
    # Create the output shapefile
    if driver.Open(output_shapefile, 1):
        driver.DeleteDataSource(output_shapefile)
    out_datasource = driver.CreateDataSource(output_shapefile)
    out_layer = out_datasource.CreateLayer("diagonal_lines", geom_type=ogr.wkbLineString, srs=source_spatial_ref)

    # Add an ID field to the output shapefile
    id_field = ogr.FieldDefn("ID", ogr.OFTInteger)
    out_layer.CreateField(id_field)

    # Loop through each feature in the input shapefile
    for feature in layer:
        geom = feature.GetGeometryRef()

        centroid = geom.Centroid()  # Get the centroid of the geometry

        # Step 1: Calculate the diagonal line length based on the size of the feature
        extent = geom.GetEnvelope()  # Bounding box (minX, maxX, minY, maxY)
        width = extent[1] - extent[0]  # XMax - XMin
        height = extent[3] - extent[2]  # YMax - YMin
        diagonal_length = math.sqrt(width ** 2 + height ** 2) * 0.5  # Half the diagonal of the bounding box

        # Step 2: Set your desired angle (in degrees) with respect to the Y-axis

        angle_rad = math.radians(angle_deg)  # Convert angle to radians

        # Step 3: Adjust the angle relative to the Y-axis and calculate offsets
        # For an angle measured counterclockwise from the positive Y-axis:
        dx = diagonal_length * math.sin(angle_rad)  # X offset for the given angle
        dy = diagonal_length * math.cos(angle_rad)  # Y offset for the given angle

        # Calculate the start and end points of the line
        x_start = centroid.GetX() - dx
        y_start = centroid.GetY() - dy
        x_end = centroid.GetX() + dx
        y_end = centroid.GetY() + dy

        # Calculate the total line length
        line_length = math.sqrt((x_end - x_start) ** 2 + (y_end - y_start) ** 2)

        # Define interval
        interval = 6  # 6 meters

        # Calculate number of intervals
        num_points = int(line_length // interval)  # Number of intervals, floor division

        # Create the diagonal line geometry with intermediate points
        diagonal_line = ogr.Geometry(ogr.wkbLineString)

        # Interpolate points along the line
        for i in range(num_points + 1):  # Include start and end points
            t = i / num_points  # Parameter from 0 to 1
            x = x_start + t * (x_end - x_start)
            y = y_start + t * (y_end - y_start)
            diagonal_line.AddPoint(x, y)

        # Step 3: Save the line to the output shapefile
        out_feature = ogr.Feature(out_layer.GetLayerDefn())
        out_feature.SetGeometry(diagonal_line)
        out_feature.SetField("ID", feature.GetFID())
        out_layer.CreateFeature(out_feature)

    # Cleanup
    datasource = None
    out_datasource = None

    print(f"Diagonal lines have been created and saved to {output_shapefile}")








    thre = 0.7
    # Paths to input TIFF file and shapefile with lines
    # tiff_file = r"E:\Mars_results_0906\grid_experiments\CTX_0_cu_10m3_11m0p4_bilinear.tif"
    # Open the raster dataset
    tiff_file = r"E:\Mars_results_0906\SFS_results\TIF\\CTX_" + str(study_area) + "_cu_10m4_11m0p4_dir_100.tif"
    sfs_depth=[]
    t

    # tiff_file = r"E:\Mars_results_0906\SFS_results\TIF\\F20_043425_2207_XN_40N354W_resample.tif"
    raster = gdal.Open(tiff_file)
    band = raster.GetRasterBand(1)  # Assume we're reading the first band

    # Get the geotransform to convert coordinates
    geotransform = raster.GetGeoTransform()
    Xsize = raster.RasterXSize
    Ysize = raster.RasterYSize


    filtered_gdf = gpd.read_file( r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  ".shp")
    filtered_gdf_line = gpd.read_file(r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_line.shp")

    to_keep = []
    for i, geom in enumerate(filtered_gdf.geometry):
        # 检查当前要素是否被其他要素包含至少 95%
        contained = False
        for j, other_geom in enumerate(filtered_gdf.geometry):
            if i != j:
                # 计算两个几何要素的交集
                intersection_area = geom.intersection(other_geom).area
                # 如果交集面积大于或等于 95% 当前要素面积，表示该要素被包含
                if intersection_area >= 0.7 * geom.area:
                    contained = True
                    break
        # 如果当前要素没有被包含，保留该要素
        to_keep.append(not contained)

    # 3. 过滤要素
    filtered_gdf = filtered_gdf[to_keep]
    filtered_gdf_line = filtered_gdf_line[to_keep]

    # 保存为新的 Shapefile
    # output_shapefile = r"E:\Mars_results_0906\process_files\buffer\buffer_result\\"+str(study_area)+"\\buffer_result_remove_"+str(thre)+"_ratio_DEM_large.shp"
    # output_shapefile = r"E:\Mars_results_0906\process_files\buffer\buffer_result\\"+str(study_area)+"\\buffer_result_remove_"+str(0.5)+".shp"
    output_shapefile = r"E:\Mars_results_0906\SFS_results\post_process\\" + str(study_area) + "\\buffer_result_" + str(
        study_area) + "_large.shp"
    output_shapefile_line = r"E:\Mars_results_0906\SFS_results\post_process\\" + str(
        study_area) + "\\buffer_result_diagonal_lines_" + str(study_area) + "_large.shp"

    filtered_gdf.to_file(output_shapefile, driver="ESRI Shapefile")
    filtered_gdf_line.to_file(output_shapefile_line, driver="ESRI Shapefile")

    # Open the shapefile containing the lines
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(output_shapefile_line, 0)  # Read-only
    layer = datasource.GetLayer()

    # Loop through each line in the shapefile
    line_values = []
    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry.GetGeometryType() == ogr.wkbLineString:
            # Extract points along the line
            points = geometry.GetPoints()  # List of (x, y) tuples
            values = [extract_pixel_value(x, y, band,geotransform,Xsize,Ysize) for x, y in points]  # Extract values for each point
            line_values.append(values)
        #print(line_values)

    # Close the datasets
    datasource = None
    raster = None

    # Print or save the extracted values
    extracted_list = []
    for i, values in tqdm(enumerate(line_values)):
        # plt.figure()
        contains_nan = np.isnan(values).any()
        if (not contains_nan):
            detrended_values = detrend(values)
            # detrended_values = values
            detrended_values = (detrended_values - np.min(detrended_values))
            # plt.plot(range(len(values)),values)
            # plt.title(str(i))
            # plt.savefig(r"E:\Mars_results_0906\process_files\buffer\lines_noDEM\18_new\\"+str(i)+".png")
            # Fit the curve
            # Fit the curve
            x_data = np.array(range(len(values)))
            x_data = (x_data - np.min(x_data))


            # print(x_data)
            # print(detrended_values)

            try:
                lower_bounds = [0, 0, 0.1, 0, np.max(x_data) / 2, 0.1]  # No constraint on b
                upper_bounds = [np.max(detrended_values), np.max(x_data) / 2, np.inf, np.max(detrended_values),
                                np.max(x_data), np.inf]

                params, params_covariance = curve_fit(fixed_function, x_data, detrended_values,
                                                      bounds=(lower_bounds, upper_bounds), maxfev=10000)
                # print(str(i))
                # print(params)

                # Calculate the fitted y values
                y_fitted = fixed_function(x_data, *params)

                # Calculate R²
                ss_res = np.sum((detrended_values - y_fitted) ** 2)
                ss_tot = np.sum((detrended_values - np.mean(detrended_values)) ** 2)
                r_squared = 1 - (ss_res / ss_tot)
                #plt.plot(x_data, y_fitted)
                a1, b1, c1, a2, b2, c2 = params
                b_center = (b1 + b2) / 2
                y_lower = fixed_function(b_center, *params)
                height = max(a1 - y_lower, a2 - y_lower)

                width = np.abs(b1 - b2)
                ratio = height * 100 / width
                b_max = max(b1, b2)
                b_min = min(b1, b2)
                #print((b_min+b_max)/2)
                #dep=((a1 + a2)/2)-y_lower
                dep = ((a1 + a2) / 2) - y_lower
                sfs_depth.append(dep)


            except:
                print("no fitting" + str(i))
    print(sfs_depth)




    # Example: 45 degrees, but you can change this value
    input_shapefile = r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_pro.shp"
    output_shapefile = r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_pro_line.shp"

    # Open the input shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(input_shapefile, 0)  # Read-only
    layer = datasource.GetLayer()
    source_spatial_ref = layer.GetSpatialRef()
    # Create the output shapefile
    if driver.Open(output_shapefile, 1):
        driver.DeleteDataSource(output_shapefile)
    out_datasource = driver.CreateDataSource(output_shapefile)
    out_layer = out_datasource.CreateLayer("diagonal_lines", geom_type=ogr.wkbLineString, srs=source_spatial_ref)

    # Add an ID field to the output shapefile
    id_field = ogr.FieldDefn("ID", ogr.OFTInteger)
    out_layer.CreateField(id_field)

    # Loop through each feature in the input shapefile
    for feature in layer:
        geom = feature.GetGeometryRef()

        centroid = geom.Centroid()  # Get the centroid of the geometry

        # Step 1: Calculate the diagonal line length based on the size of the feature
        extent = geom.GetEnvelope()  # Bounding box (minX, maxX, minY, maxY)
        width = extent[1] - extent[0]  # XMax - XMin
        height = extent[3] - extent[2]  # YMax - YMin
        diagonal_length = math.sqrt(width ** 2 + height ** 2) * 0.5  # Half the diagonal of the bounding box

        # Step 2: Set your desired angle (in degrees) with respect to the Y-axis

        angle_rad = math.radians(angle_deg)  # Convert angle to radians

        # Step 3: Adjust the angle relative to the Y-axis and calculate offsets
        # For an angle measured counterclockwise from the positive Y-axis:
        dx = diagonal_length * math.sin(angle_rad)  # X offset for the given angle
        dy = diagonal_length * math.cos(angle_rad)  # Y offset for the given angle

        # Calculate the start and end points of the line
        x_start = centroid.GetX() - dx
        y_start = centroid.GetY() - dy
        x_end = centroid.GetX() + dx
        y_end = centroid.GetY() + dy

        # Calculate the total line length
        line_length = math.sqrt((x_end - x_start) ** 2 + (y_end - y_start) ** 2)

        # Define interval
        interval = 6  # 6 meters

        # Calculate number of intervals
        num_points = int(line_length // interval)  # Number of intervals, floor division

        # Create the diagonal line geometry with intermediate points
        diagonal_line = ogr.Geometry(ogr.wkbLineString)

        # Interpolate points along the line
        for i in range(num_points + 1):  # Include start and end points
            t = i / num_points  # Parameter from 0 to 1
            x = x_start + t * (x_end - x_start)
            y = y_start + t * (y_end - y_start)
            diagonal_line.AddPoint(x, y)

        # Step 3: Save the line to the output shapefile
        out_feature = ogr.Feature(out_layer.GetLayerDefn())
        out_feature.SetGeometry(diagonal_line)
        out_feature.SetField("ID", feature.GetFID())
        out_layer.CreateFeature(out_feature)

    # Cleanup
    datasource = None
    out_datasource = None

    print(f"Diagonal lines have been created and saved to {output_shapefile}")












    sfs_depth_2=[]
    tiff_file =r"E:\Mars_results_0906\DEM_data\CTX_"+str(study_area)+"_hirise1.tif"

    # tiff_file = r"E:\Mars_results_0906\SFS_results\TIF\\F20_043425_2207_XN_40N354W_resample.tif"
    raster = gdal.Open(tiff_file)
    band = raster.GetRasterBand(1)  # Assume we're reading the first band



    # Get the geotransform to convert coordinates
    geotransform = raster.GetGeoTransform()
    Xsize = raster.RasterXSize
    Ysize = raster.RasterYSize


    filtered_gdf = gpd.read_file( r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_pro.shp")
    filtered_gdf_line = gpd.read_file(r"E:\Mars_results_0906\DEM_data\clip_" + str(study_area) +  "_pro_line.shp")

    to_keep = []
    for i, geom in enumerate(filtered_gdf.geometry):
        # 检查当前要素是否被其他要素包含至少 95%
        contained = False
        for j, other_geom in enumerate(filtered_gdf.geometry):
            if i != j:
                # 计算两个几何要素的交集
                intersection_area = geom.intersection(other_geom).area
                # 如果交集面积大于或等于 95% 当前要素面积，表示该要素被包含
                if intersection_area >= 0.7 * geom.area:
                    contained = True
                    break
        # 如果当前要素没有被包含，保留该要素
        to_keep.append(not contained)

    # 3. 过滤要素
    filtered_gdf = filtered_gdf[to_keep]
    filtered_gdf_line = filtered_gdf_line[to_keep]

    # 保存为新的 Shapefile
    # output_shapefile = r"E:\Mars_results_0906\process_files\buffer\buffer_result\\"+str(study_area)+"\\buffer_result_remove_"+str(thre)+"_ratio_DEM_large.shp"
    # output_shapefile = r"E:\Mars_results_0906\process_files\buffer\buffer_result\\"+str(study_area)+"\\buffer_result_remove_"+str(0.5)+".shp"
    output_shapefile = r"E:\Mars_results_0906\SFS_results\post_process\\" + str(study_area) + "\\buffer_result_" + str(
        study_area) + "_large.shp"
    output_shapefile_line = r"E:\Mars_results_0906\SFS_results\post_process\\" + str(
        study_area) + "\\buffer_result_diagonal_lines_" + str(study_area) + "_large.shp"

    filtered_gdf.to_file(output_shapefile, driver="ESRI Shapefile")
    filtered_gdf_line.to_file(output_shapefile_line, driver="ESRI Shapefile")

    # Open the shapefile containing the lines
    driver = ogr.GetDriverByName("ESRI Shapefile")
    datasource = driver.Open(output_shapefile_line, 0)  # Read-only
    layer = datasource.GetLayer()

    # Loop through each line in the shapefile
    line_values = []
    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry.GetGeometryType() == ogr.wkbLineString:
            # Extract points along the line
            points = geometry.GetPoints()  # List of (x, y) tuples
            values = [extract_pixel_value(x, y, band,geotransform,Xsize,Ysize) for x, y in points]  # Extract values for each point
            line_values.append(values)

    # Close the datasets
    datasource = None
    raster = None

    # Print or save the extracted values
    extracted_list = []
    for i, values in tqdm(enumerate(line_values)):
        # plt.figure()
        contains_nan = np.isnan(values).any()
        if (not contains_nan):
            #detrended_values = detrend(values)
            detrended_values = values
            detrended_values = (detrended_values - np.min(detrended_values))
            # plt.plot(range(len(values)),values)
            # plt.title(str(i))
            # plt.savefig(r"E:\Mars_results_0906\process_files\buffer\lines_noDEM\18_new\\"+str(i)+".png")
            # Fit the curve
            # Fit the curve
            x_data = np.array(range(len(values)))
            x_data = (x_data - np.min(x_data))


            # print(x_data)
            # print(detrended_values)

            try:
                lower_bounds = [0, 0, 0.1, 0, np.max(x_data) / 2, 0.1]  # No constraint on b
                upper_bounds = [np.max(detrended_values), np.max(x_data) / 2, np.inf, np.max(detrended_values),
                                np.max(x_data), np.inf]

                params, params_covariance = curve_fit(fixed_function, x_data, detrended_values,
                                                      bounds=(lower_bounds, upper_bounds), maxfev=10000)
                # print(str(i))
                # print(params)

                # Calculate the fitted y values
                y_fitted = fixed_function(x_data, *params)
                # plt.plot(x_data, detrended_values)
                # plt.show()


                # Calculate R²
                ss_res = np.sum((detrended_values - y_fitted) ** 2)
                ss_tot = np.sum((detrended_values - np.mean(detrended_values)) ** 2)
                r_squared = 1 - (ss_res / ss_tot)
                #plt.plot(x_data, y_fitted)
                a1, b1, c1, a2, b2, c2 = params
                b_center = (b1 + b2) / 2
                y_lower = fixed_function(b_center, *params)
                height = max(a1 - y_lower, a2 - y_lower)

                width = np.abs(b1 - b2)
                ratio = height * 100 / width
                b_max = max(b1, b2)
                b_min = min(b1, b2)
                print(a1,a2)
                #print((a1 - a2)/2-y_lower)
                #print(a1,a2,y_lower)
                if (max(detrended_values)-min(detrended_values))>2000:
                    sfs_depth_2.append(np.nan)
                else:
                    sfs_depth_2.append(max(detrended_values)-min(detrended_values))

            except:
                print("no fitting" + str(i))
                sfs_depth_2.append(np.nan)
    #print(sfs_depth_2)
    #print(sfs_depth)
    #print(np.array(sfs_depth)-np.array(sfs_depth_2))
    diff=abs(np.array(sfs_depth)-np.array(sfs_depth_2))
    sfs_depth_2=np.array(sfs_depth_2)
    sfs_depth_filter=sfs_depth_2[diff>-2000]
    print(diff[diff>-2000])
    diff=diff[diff>-2000]
    #diff = diff[diff<20]
    sfs_depth=np.array(sfs_depth)
    sfs_depth_2 = np.array(sfs_depth_2)


    #print(abs(diff))
    print(np.average(diff))
    print(sfs_depth_filter)
    print(np.max(diff))
    #rmse = np.sqrt(np.mean((sfs_depth[sfs_depth_2>-2000] - sfs_depth_2[sfs_depth_2>-2000]) ** 2))
    print(len(diff[diff<5])/len(diff))
    #print(np.average(diff/sfs_depth[sfs_depth>0]))
    plt.hist(diff,bins=50)
    plt.show()