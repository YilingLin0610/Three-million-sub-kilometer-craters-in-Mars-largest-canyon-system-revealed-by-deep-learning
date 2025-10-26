from merge import post_process
import os
import time
def list_dirs(path="."):
    """只列出目录，不列文件"""
    return [
        name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))
    ]
if __name__ == "__main__":
    loop=1
    while(loop==1):
        total_paths = [
                       r"/data/B/Yiling/Mars_SFS/CTX_images_process_1/",
                       r"/data/B/Yiling/Mars_SFS/CTX_images_process/",
                       r"/data/B/Yiling/Mars_SFS/CTX_images_process_2/",
                       r"/data/B/Yiling/Mars_SFS/CTX_images_process_3/",
                       r"/data/B/Yiling/Mars_SFS/CTX_images_process_4/",
            r"/data/C/Yiling/Mars_SFS/CTX_images_process_1/",
                       ]
        for total_path in total_paths:
            print(list_dirs(total_path))
            for mosaic_ID in list_dirs(total_path):
                if (os.path.exists(total_path + r"/" + mosaic_ID + "/result_data")):
                    result_data_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/result_data"))
                    test_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/test"))
                    print(total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83.shp")
                    if (os.path.exists(total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83.shp")):
                        print(mosaic_ID + "processed")
                    elif (
                            os.path.exists(
                                total_path + r"/" + mosaic_ID + "/merged_result_update_0.5_roundness_0.83_pro.shp")):
                        print(mosaic_ID + "processed")
                    elif (result_data_length == test_length):
                        print(mosaic_ID + "processing")
                        try:
                            SFS_path = total_path + mosaic_ID + "/"
                            post_process(mosaic_ID, SFS_path, 0.5, 0.83)
                        except Exception as e:
                            # optional last-resort guardrail
                            print(f"Unexpected: {e!r}")

    time.sleep(600)
    print("sleeping")
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
