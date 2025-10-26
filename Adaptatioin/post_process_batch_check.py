from merge_tif import post_process
import os
def list_dirs(path="."):
    """只列出目录，不列文件"""
    return [
        name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))
    ]
if __name__ == "__main__":
    total_paths = [r"G:\Mars\Mars_SFS_batch\CTX_images_process_C_1\\"]
    for total_path in total_paths:
        print(list_dirs(total_path))
        for mosaic_ID in ["MurrayLab_GlobalCTXMosaic_V01_E100_N-68"]:
            result_data_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/result_data"))
            #test_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/test"))
            print(total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83.shp")
            if (os.path.exists(total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83k.shp")):
                print(mosaic_ID + "processed")
            # elif (os.path.exists(total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83_pro.shp")):
            #     print(mosaic_ID + "processed")
            else:
                print(mosaic_ID + "processing")
                SFS_path = total_path + mosaic_ID + "/"
                post_process(mosaic_ID, SFS_path, 0.5, 0.83)
    #total_path=r"/data/C/Yiling/Mars_SFS/CTX_images_process_1/"
    #os.listdir(total_path)


            #except Exception as e:
            #    print(e)
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E148_N-88'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path,0.5,0.83)
    #
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-048_N-08'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path, 0.5, 0.7)
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N00'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path, 0.5, 0.7)
    # post_process(mosaic_ID, SFS_path, 0.5, 0.75)
    # post_process(mosaic_ID, SFS_path,0.5,0.85)
    # post_process(mosaic_ID, SFS_path, 0.5, 0.9)


    # MurrayLab_CTX_V01_E104_N - 04
    # _SeamMap_merged.shp
