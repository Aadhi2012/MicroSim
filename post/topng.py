import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

def has_h5_extension(filename):
    _, file_extension = os.path.splitext(filename)
    return file_extension.lower() == '.h5'

def list_h5_files(directory):
    h5_files_list = []
    for item in os.listdir(directory):
        full_path = os.path.join(directory, item)
        if os.path.isfile(full_path):
            if(has_h5_extension(full_path)):
                h5_files_list.append(full_path)
    return h5_files_list

def main():
    if len(sys.argv) > 1 :
        main_dir = sys.argv[1]
    else:
        main_dir = "DATA"

    files_h5 = list_h5_files(main_dir)

    if(len(sys.argv) > 2 ):
        plotKey = sys.argv[2]
    else:
        plotKey = "all"
    
    if(len(sys.argv) > 3 ):
        timeStep = sys.argv[3]
    else:
        timeStep = "all"
    
    for filepath in files_h5:
        if  ("_" + timeStep + "." in filepath) or timeStep == "all" :
            with h5py.File(filepath, "r") as DATA:
                KEYS = list(DATA.keys())
                print("\n- File %s has keys: %s" % (filepath , KEYS))
                for key in KEYS:
                    if key == plotKey or plotKey == "all" :
                        #print(f"- {key} found in file {filepath}")
                        npDATA = DATA[key][()]
                        npShape = np.shape(npDATA)
                        plt.figure()
                        print(f"- shape of dataset : {npShape}")
                        plt.imshow( npDATA) #[int(npShape[0]/2)] )
                        plt.xlabel("Y-axis")
                        plt.ylabel("Z-axis")
                        plt.colorbar()
                        if (len(sys.argv) > 4) and sys.argv[4] == "show":
                            plt.show()
                        else:
                            savename = filepath.replace(".h5","_") + key + ".png"
                            plt.savefig(savename)
                            print(f"- {savename} saved!")
                        plt.clf()
                        plt.close()
                    else:
                        continue
        else:
            continue


if __name__ == "__main__":
    main()
