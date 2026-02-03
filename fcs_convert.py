import sys
import fcsparser

if __name__ == "__main__":
    try:
        fp_load = sys.argv[1]
        print("reading fcs file: ", fp_load)
        meta, data = fcsparser.parse(fp_load)
        fp_save = fp_load.split(".")[0]+".csv"
        print("save csv file: ",fp_save) 
        data.to_csv(fp_save)
    except Exception as err:
        print("ERROR: ", err)

