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


def write_xmf(fname, npShape, keys):
    f = open(fname,'w')
    fname = fname.replace(".xmf", ".h5")
    shape_str = ""
    print(np.shape(npShape))
    if np.shape(npShape)[0] == 2 :
    	shape_str = str(npShape[0]) + " " + str(npShape[1])
    else:
    	shape_str = str(npShape[0]) + " " + str(npShape[1]) + " " + str(npShape[2])
    f.write("<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf Version=\"2.0\">\n<Domain>\n<Grid Name=\"mesh1\" GridType=\"Uniform\">")
    f.write("<Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"" + shape_str + "\"/>")
    f.write("<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n<DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n 0 0 0 \n</DataItem>\n<DataItem Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n 1.000000 1.000000 1.000000\n</DataItem>\n</Geometry>")
    f.write("\n\n")

    for key in keys:
        f.write("<Attribute Name=\" "+key+"\" AttributeType=\"Scalar\" Center=\"Node\">\n<DataItem Dimensions=\""+shape_str+"\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n")
        f.write(fname + ":/" + key)
        f.write("\n</DataItem>\n </Attribute>\n")
        f.write("\n")

    f.write("\n")
    f.write("</Grid>\n</Domain>\n</Xdmf>")
    f.close()



def main():
    if len(sys.argv) > 1 :
        main_dir = sys.argv[1]
    else:
        main_dir = "DATA"

    files_h5 = list_h5_files(main_dir)

    if(len(sys.argv) > 2 ):
        timeStep = sys.argv[2]
    else:
        timeStep = "all"

    for filepath in files_h5:
        xmffilename = filepath.replace(".h5", ".xmf")
        
        if os.path.isfile(xmffilename):
            print(f"- file {xmffilename} exists!")
            continue

        if  timeStep  in filepath or timeStep == "all" :
            with h5py.File(filepath, "r") as DATA:
                KEYS = list(DATA.keys())
                print("\n- File %s has keys: %s" % (filepath , KEYS))
                for key in KEYS:
                    npShape = np.shape(DATA[key][()])
                
                write_xmf(xmffilename, npShape, KEYS)
                DATA.close()

        else:
            continue
    

if __name__ == "__main__":
    main()
