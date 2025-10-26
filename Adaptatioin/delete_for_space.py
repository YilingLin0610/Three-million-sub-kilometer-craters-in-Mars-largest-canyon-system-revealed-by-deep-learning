

from merge import post_process
import os
from post_process_batch import list_dirs
import shutil
from tqdm import tqdm
if __name__ == "__main__":
    # total_paths = [r"/data/B/Yiling/Mars_SFS/CTX_images_process/",
    #                r"/data/B/Yiling/Mars_SFS/CTX_images_process_1/",
    #                r"/data/B/Yiling/Mars_SFS/CTX_images_process_2/",
    #                r"/data/B/Yiling/Mars_SFS/CTX_images_process_3/",
    #                r"/data/B/Yiling/Mars_SFS/CTX_images_process_4/",r"/data/B/Yiling/Mars_SFS/CTX_images_process_5/",
    #              r"/data/B/Yiling/Mars_SFS/CTX_images_process_6/",r"/data/C/Yiling/Mars_SFS/CTX_images_process_1/",
    #              r"/data/C/Yiling/Mars_SFS/CTX_images_process_2/",r"/data/C/Yiling/Mars_SFS/CTX_images_process_3/",
    #              r"/data/C/Yiling/Mars_SFS/CTX_images_process_4/",r"/data/C/Yiling/Mars_SFS/CTX_images_process_5/",
    #              r"/data/C/Yiling/Mars_SFS/CTX_images_process_6/"]
    total_paths = [r"/data/B/Yiling/Mars_SFS/CTX_images_process_6/"]
    for total_path in total_paths:
        print(list_dirs(total_path))
        for mosaic_ID in tqdm(list_dirs(total_path)):


            # if (os.path.exists(total_path + r"/" + mosaic_ID + "/result_data") and os.path.exists(total_path + r"/" + mosaic_ID + "/test")):
            #     result_data_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/result_data"))
            #     test_length = len(os.listdir(total_path + r"/" + mosaic_ID + "/test"))
            #     #if(result_data_length==test_length):
            #     if (result_data_length == test_length):
            #         print(mosaic_ID + "processing")
            #
            #         SFS_path = total_path + mosaic_ID + "/"
            #         path1 = SFS_path + "tifs/"
            #         shutil.rmtree(path1)
            #         path2 = SFS_path + "tifs_jpgs/"
            #         shutil.rmtree(path2)
            print(total_path + r"/" + mosaic_ID + "/results/merge_prediction_shp"+mosaic_ID+".shp")
            if (os.path.exists(total_path + r"/" + mosaic_ID + "/results/merge_prediction_shp//"+mosaic_ID+".shp")):

                print(mosaic_ID)
                SFS_path = total_path + mosaic_ID + "/"
                try:
                    path1 = SFS_path + "tifs/"
                    shutil.rmtree(path1)
                except:
                    pass
                try:
                    path2 = SFS_path + "tifs_jpgs/"
                    shutil.rmtree(path2)
                except:
                    pass
                try:
                    path3 = SFS_path + "results/merge_prediction_tifs"
                    print(path3)
                    shutil.rmtree(path3)
                except:
                    pass
            # if (os.path.exists(total_path + r"/" + mosaic_ID + "/result_data")):
            #     print(mosaic_ID)
            #     SFS_path = total_path + mosaic_ID + "/"
            #     try:
            #         path1 = SFS_path +  "/result_data"
            #         shutil.rmtree(path1)
            #
            #     except:
            #         pass


