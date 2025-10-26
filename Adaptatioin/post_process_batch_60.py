from merge_60 import post_process
import os
from tqdm import tqdm
print("yes")
if __name__ == "__main__":
    total_path=r"/data/A/Yiling/Mars_SAM_crater/CTX_images_Valles//"
    Mosaic_names = [ f for f in os.listdir(total_path) if os.path.isdir(os.path.join(total_path, f))]
    #Mosaic_names=['MurrayLab_GlobalCTXMosaic_V01_E-052_N00']
    for mosaic_ID in tqdm(Mosaic_names):
        if(os.path.exists( total_path + mosaic_ID + "/"+"merged_result_update_60_0.85.shp")):
            print(mosaic_ID+"processed")
        else:
            #try:
            SFS_path = total_path + mosaic_ID + "/"
            print(SFS_path)
            post_process(mosaic_ID, SFS_path)
            # except:
            #     print(mosaic_ID+"failed")




    #mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N-04'

    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-036_N-20'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)
    # mosaic_ID='MurrayLab_GlobalCTXMosaic_V01_E-088_N00'
    # SFS_path = total_path + mosaic_ID + "/"
    # post_process(mosaic_ID, SFS_path)


    # MurrayLab_CTX_V01_E104_N - 04
    # _SeamMap_merged.shp
