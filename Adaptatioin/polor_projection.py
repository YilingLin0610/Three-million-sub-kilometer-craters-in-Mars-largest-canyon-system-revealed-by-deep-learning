
from osgeo import gdal, ogr
from typing import Optional, Tuple

def reproject_raster_to_polar(

    in_tif: str,
    out_tif: str,
    pole,                 # "north" 或 "south"
    lon0: float = 0.0,                   # 极投影中央经线（度）
    src_srs: Optional[str] = "+proj=eqc +lat_ts=0 +lon_0=0 +R=3396190 +units=m +no_defs",  # 源投影（None=读文件）；建议用 PROJ4
    lat_ts: Optional[float] = None,      # Variant B：标准纬线；None 则用 k（Variant A）
    k: float = 1.0,                      # Variant A：比例因子
    R: float = 3396190.0,                # 火星球半径（米）
    target_res: Optional[Tuple[float, float]] = None,  # (xRes, yRes)，不传则让 GDAL 自动
    resample: str = "bilinear",          # 连续数据用 bilinear/cubic；分类用 near
    src_nodata: Optional[float] = None,
    dst_nodata: Optional[float] = None,
    # 若提供 ref_tif：继承其 投影+范围+像元数（完全对齐到参考格网）
    ref_tif: Optional[str] = None
) -> str:
    """
    将 in_tif 重投影为极地立体（Polar Stereographic）。
    若传入 ref_tif，则输出将与 ref_tif 完全对齐（投影、范围、分辨率、行列数）。
    返回 out_tif 路径。
    """
    # 读输入

    src_ds = gdal.Open(in_tif)
    if src_ds is None:
        raise RuntimeError(f"无法打开输入文件：{in_tif}")

    # 源投影：未显式给时，尝试读文件
    if src_srs is None:
        src_proj = src_ds.GetProjection()
        if not src_proj or not src_proj.strip():
            raise RuntimeError("输入无投影定义，请传入 src_srs（建议用 PROJ4）。")
        src_srs = src_proj

    # 目标投影：极地立体（优先 Variant B，其次 Variant A）
    if pole > 0:
        lat0 = 90
    else:
        lat0 = -90
    if (lat_ts is not None) and (abs(lat_ts) >= 90):
        raise ValueError("lat_ts 不能等于 ±90；请给如 70 或 -71 这样的值。")

    if lat_ts is not None:
        dst_srs = f"+proj=stere +lat_0={lat0} +lat_ts={lat_ts} +lon_0={lon0} +R={R} +units=m +no_defs"
    else:
        dst_srs = f"+proj=stere +lat_0={lat0} +lon_0={lon0} +k={k} +R={R} +units=m +no_defs"

    # 组装 Warp 选项
    warp_kwargs = dict(
        srcSRS=src_srs,
        dstSRS=dst_srs,
        resampleAlg=resample,
        format="GTiff",
        creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]
    )
    if src_nodata is not None:
        warp_kwargs["srcNodata"] = src_nodata
    if dst_nodata is not None:
        warp_kwargs["dstNodata"] = dst_nodata
    if target_res is not None and ref_tif is None:
        warp_kwargs["xRes"], warp_kwargs["yRes"] = target_res

    # 若提供参考栅格：完全对齐到参考格网（并继承其投影！）
    if ref_tif is not None:
        b = gdal.Open(ref_tif)
        if b is None:
            raise RuntimeError(f"无法打开参考文件：{ref_tif}")
        gt = b.GetGeoTransform()
        width, height = b.RasterXSize, b.RasterYSize
        xmin = gt[0]
        xmax = xmin + gt[1] * width
        ymax = gt[3]
        ymin = ymax + gt[5] * height

        # 注意：这里选择**继承参考栅格的投影**与格网
        # 若你的 ref_tif 就是极地立体，这样可保证 1:1 对齐。
        warp_kwargs.update(
            dstSRS=b.GetProjection() or dst_srs,  # 若 ref 无投影，则退回到构造的极立体
            outputBounds=[xmin, ymin, xmax, ymax],
            width=width,
            height=height
        )

    # 执行重投影
    print("源投影:", src_srs)
    print("目标投影:", warp_kwargs.get("dstSRS"))
    gdal.Warp(out_tif, in_tif, **warp_kwargs)
    print("完成:", out_tif)
    return out_tif

# ======= 必改：输入/输出路径 =======
IN_SHP  = r"G:\Mars\Mars_SFS_batch\CTX_images_process_C_1\MurrayLab_GlobalCTXMosaic_V01_E076_N76\MurrayLab_GlobalCTXMosaic_V01_E076_N76_buffer_split.shp"     # 原始 Shapefile
OUT_SHP = r"G:\Mars\Mars_SFS_batch\CTX_images_process_C_1\MurrayLab_GlobalCTXMosaic_V01_E076_N76\MurrayLab_GlobalCTXMosaic_V01_E076_N76_buffer_split_pro.shp"  # 输出 Shapefile

# ======= 选项：北极 or 南极、中央经线 =======
POLE = "south"   # "north" 或 "south"
LON0 = 0.0       # 极投影中央经线（度），按你的工程需要设（如 0、90、180、16.281 等）
from osgeo import gdal
from typing import Optional  # 3.8 用 Optional/Union
# 常用：火星等距圆柱（Equirectangular, clon=0, ocentric, IAU2015 球半径）
MARS_EQC_CLON0 = "+proj=eqc +lat_ts=0 +lon_0=0 +R=3396190 +units=m +no_defs"

def reproject_shp_to_polar(
    in_shp: str,
    out_shp: str,
    pole,                 # "north" 或 "south"
    lon0: float = 0.0,                   # 极投影中央经线（度）
    src_srs: Optional[str] = MARS_EQC_CLON0,  # 源投影；None 表示从 .prj 读取
    lat_ts: Optional[float] = None,      # 标准纬线（度）；None 则用 k（Variant A）
    k: float = 1.0,                      # 比例因子（Variant A）
    R: float = 3396190.0,                # 火星球半径（米）
    make_valid: bool = True,
    skip_failures: bool = True,
    out_format: str = "ESRI Shapefile",
    encoding: str = "UTF-8"
) -> str:
    """将 Shapefile 重投影为火星极地立体投影（Polar Stereographic）。"""
    # 目标极：北极或南极
    if pole>0:
        lat0 = 90
    else:
        lat0 = -90

    #lat0 = 90 if str(pole).lower().startswith("n") else -90

    # 目标投影（Variant B 优先，其次 Variant A）
    if lat_ts is not None:
        dst_srs = f"+proj=stere +lat_0={lat0} +lat_ts={lat_ts} +lon_0={lon0} +R={R} +units=m +no_defs"
    else:
        dst_srs = f"+proj=stere +lat_0={lat0} +lon_0={lon0} +k={k} +R={R} +units=m +no_defs"

    # 若未显式提供源投影，则尝试从 .prj 读取
    if src_srs is None:
        ds_in = gdal.OpenEx(in_shp)
        if ds_in is None:
            raise RuntimeError(f"无法打开输入：{in_shp}")
        lyr = ds_in.GetLayer(0)
        sref = lyr.GetSpatialRef() if lyr is not None else None
        wkt = sref.ExportToWkt() if sref is not None else None
        if not wkt or not wkt.strip():
            raise RuntimeError("输入 Shapefile 无投影定义且未设置 src_srs，请显式指定源投影。")
        src_srs = wkt

    opts = gdal.VectorTranslateOptions(
        dstSRS=dst_srs,
        srcSRS=src_srs,
        format=out_format,
        layerCreationOptions=[f"ENCODING={encoding}"],
        makeValid=make_valid,
        skipFailures=skip_failures
    )

    gdal.VectorTranslate(destNameOrDestDS=out_shp, srcDS=in_shp, options=opts)
    return out_shp





