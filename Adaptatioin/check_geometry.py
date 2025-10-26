import re
import requests
from bs4 import BeautifulSoup
import os
import geopandas as gpd
from tqdm import tqdm
import pandas as pd
import time
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






def fetch_ctx_metadata(url: str, retries: int = 3, backoff: float = 1.0, timeout: int = 10):
    """
    最多重试 retries 次；指数退避：backoff, 2*backoff, 4*backoff ...
    失败最终抛异常，由外层记录后继续。
    """
    headers = {"User-Agent": "Mozilla/5.0"}
    last_err = None
    for i in range(1, retries + 1):
        try:
            resp = requests.get(url, headers=headers, timeout=timeout)
            resp.raise_for_status()
            soup = BeautifulSoup(resp.text, "html.parser")
            text = soup.get_text(" ", strip=True)

            fields = ["Subsolar Longitude", "Subsolar Latitude", "Center Lat", "Center Lon"]
            results = {}
            for field in fields:
                m = re.search(rf"{field}\s*[:=]?\s*([-+]?\d+(?:\.\d+)?)", text)
                results[field] = float(m.group(1)) if m else None
            return results
        except requests.RequestException as e:
            last_err = e
            if i < retries:
                time.sleep(backoff * (2 ** (i - 1)))  # 1x, 2x, 4x ...
            else:
                # 最后一次失败，抛给外层
                raise last_err

def batch_check(total_path):
    records = []
    failures = []   # ← 新增：记录失败

    for mosaic_ID in tqdm(list_dirs(total_path)):
        SFS_path = total_path + mosaic_ID + "/"
        splits = mosaic_ID.split("_")
        latitude = int(splits[-1][1:])

        merge_shapefile_path = (
            SFS_path + splits[0] + "_CTX_" + splits[2] + "_" + splits[3] + "_" + splits[4] + "_SeamMap_merged.shp"
        )
        gdf = gpd.read_file(merge_shapefile_path)

        for index, row in gdf.iterrows():
            name = row.PRODUCT_ID
            url = f"https://viewer.mars.asu.edu/viewer/ctx/{name}#T=2&P={name}"
            try:
                # 带重试的抓取（可按需改重试参数）
                results = fetch_ctx_metadata(url, retries=3, backoff=1.0, timeout=10)
            except Exception as e:
                # 失败就记录，并继续下一个
                failures.append({
                    "Head Folder": mosaic_ID,
                    "Sub Folder": name,
                    "URL": url,
                    "Error": str(e),
                })
                continue

            results_match = compare_and_export_to_excel(row, results, mosaic_ID, name, url)
            if results_match is not None:
                print(mosaic_ID)
                records.append(results_match)

    # 原有导出
    df = pd.DataFrame(records)
    df.to_excel(total_path + "geometry_check.xlsx", index=False)

    # 失败清单导出（若有）
    if failures:
        pd.DataFrame(failures).to_excel(total_path + "geometry_fetch_failures.xlsx", index=False)
        # 或者 CSV：
        # pd.DataFrame(failures).to_csv(total_path + "geometry_fetch_failures.csv", index=False)



if __name__ == "__main__":

    total_paths=[

                 r"/data/C/Yiling/Mars_SFS/CTX_images_process_6/"]
    for total_path in total_paths:
        try:
            batch_check(total_path)
        except:
            pass














#
# url = "https://viewer.mars.asu.edu/viewer/ctx/K11_057729_1101_XN_69S236W#T=2&P=K11_057729_1101_XN_69S236W"  # 替换成实际网页地址
#
# # 请求网页
# headers = {"User-Agent": "Mozilla/5.0"}
# resp = requests.get(url, headers=headers)
# resp.raise_for_status()
#
# # 解析网页
# soup = BeautifulSoup(resp.text, "html.parser")
# text = soup.get_text(" ", strip=True)  # 把所有文字拉出来
#
# # 用正则匹配 "Subsolar Longitude" 后面的数字
# match = re.search(r"Subsolar Longitude\s*[:=]?\s*([-+]?\d+(?:\.\d+)?)", text)
#
# if match:
#     value = float(match.group(1))
#     print("Subsolar Longitude =", value)
# else:
#     print("没有找到 Subsolar Longitude")







