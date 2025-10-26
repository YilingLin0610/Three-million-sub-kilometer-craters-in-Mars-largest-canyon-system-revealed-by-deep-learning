import re
import requests
from bs4 import BeautifulSoup
import os
import geopandas as gpd
from tqdm import tqdm
import pandas as pd
from pathlib import Path
import math
import numpy as np
# 映射关系：网页字段名 -> row中的字段名
FIELD_MAP = {
    "Center Lon": "CLONG",
    "Center Lat": "CLAT",
    "Subsolar Longitude": "SB_SLR_LNG",
    "Subsolar Latitude": "SB_SLR_LAT",
}
def list_dirs(path="."):
    """只列出目录，不列文件"""
    return [
        name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))
    ]
def compare_and_export_to_excel(
    row,
    extracted: dict,
    head_folder: str,
    sub_folder: str,
    url: str,
):
    """
    四个字段逐一对比（阈值±0.15）。
    - 若全部匹配 → 返回 None
    - 若存在任意不匹配 → 返回包含三列(Head/Sub/URL) + 八个值(Web/Local×4)的一行dict
    """
    record = {
        "Head Folder": head_folder,
        "Sub Folder": sub_folder,
        "URL": url,
    }

    has_mismatch = False
    value_cols = {}

    for web_field, row_attr in FIELD_MAP.items():
        local_val = getattr(row, row_attr)
        web_val = extracted.get(web_field, None)

        prefix = web_field.replace(" ", "_")
        value_cols[f"{prefix}_Web"] = web_val
        value_cols[f"{prefix}_Local"] = local_val

        match = (web_val is not None) and (abs(web_val - local_val) <= 0.15)
        if not match:
            has_mismatch = True

    # 调试：需要的话可以打印
    # print(value_cols)

    if has_mismatch:
        record.update(value_cols)   # 合并八个值
        return record               # 返回包含 Head/Sub/URL + 八个值

    return None                     # 全部匹配返回空






def fetch_ctx_metadata(url: str):
    headers = {"User-Agent": "Mozilla/5.0"}
    resp = requests.get(url, headers=headers)
    resp.raise_for_status()

    # 解析网页
    soup = BeautifulSoup(resp.text, "html.parser")
    text = soup.get_text(" ", strip=True)

    # 定义要查找的字段
    fields = ["Subsolar Longitude", "Subsolar Latitude", "Center Lat", "Center Lon"]
    results = {}

    for field in fields:
        # 匹配数字（可以是正负、小数）
        match = re.search(rf"{field}\s*[:=]?\s*([-+]?\d+(?:\.\d+)?)", text)
        if match:
            results[field] = float(match.group(1))
        else:
            results[field] = None

    return results

def batch_check(total_path,mosaic_IDs):
    records = []
    for mosaic_ID in mosaic_IDs:
        # print(mosaic_ID)
        SFS_path = total_path + mosaic_ID + "/"
        splits = mosaic_ID.split("_")
        latitude = int(splits[-1][1:])

        merge_shapefile_path = SFS_path + splits[0] + "_CTX_" + splits[2] + "_" + splits[3] + "_" + splits[
            4] + "_SeamMap_merged.shp"
        gdf = gpd.read_file(merge_shapefile_path)
        crs = gdf.crs
        merged_gdf = gpd.GeoDataFrame()

        for idx, row in gdf.iterrows():
            name = row.PRODUCT_ID
            print(name)
            # center_lon=row.CLONG
            # center_lat=row.CLAT
            # solar_lon=row.SB_SLR_LNG
            # solar_lat=row.SB_SLR_LAT
            # print(name)
            url = "https://viewer.mars.asu.edu/viewer/ctx/" + name + "#T=2&P=" + name
            results = fetch_ctx_metadata(url)
            # 逐个字段写回（存在才改）
            v = results.get("Center Lon")
            if v is not None: gdf.at[idx, "CLONG"] = float(v)

            v = results.get("Center Lat")
            if v is not None: gdf.at[idx, "CLAT"] = float(v)

            v = results.get("Subsolar Longitude")
            if v is not None: gdf.at[idx, "SB_SLR_LNG"] = float(v)

            v = results.get("Subsolar Latitude")
            if v is not None: gdf.at[idx, "SB_SLR_LAT"] = float(v)
        out_path = Path(merge_shapefile_path).with_name(Path(merge_shapefile_path).stem + "_webfix.shp")
        gdf.to_file(out_path, driver="ESRI Shapefile", encoding="utf-8")
        print(f"✅ 已写出: {out_path}")
            #print(results)
            # # print(results)
            # results_match = compare_and_export_to_excel(row, results, mosaic_ID, name, url)
            # if results_match is None:
            #     pass
            # else:
            #     # print(results_match)
            #     #print(mosaic_ID)
            #     records.append(results_match)
    #df = pd.DataFrame(records)
    #df.to_excel(total_path + "geometry_check.xlsx")
    # —— 一次性保存 —— #




if __name__ == "__main__":

    total_paths=[
                 r"/data/C/Yiling/Mars_SFS/CTX_images_process_6/"]
    for total_path in total_paths:
        #batch_check(total_path)
        excel = total_path+r"geometry_check.xlsx"  # 改成你的路径
        df = pd.read_excel(excel)

        # 列名清理（防止有前后空格）
        df.columns = df.columns.str.strip()

        col = "Head Folder"
        uniq = (df[col]
                .dropna()
                .astype(str)
                .str.strip()
                .drop_duplicates()
                .sort_values()
                .tolist())
        print(uniq)

        print(f"共 {len(uniq)} 个唯一 Head Folder：")
        batch_check(total_path, uniq)







