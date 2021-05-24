import sys
from sequence import *

def main():

    if len(sys.argv)!=3:
        print("ERROR: Invalid number of arguments for main.py\nWrite: python main.py <input_file.dat> <output_file.dat>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

        data = np.loadtxt(input_file)
        seq = Sequence(data) 
        sort_seq = seq.sorted(key='best')
        np.savetxt(output_file, sort_seq.data[sort_seq.contour])
        print("INFO: Successfully saved sorted data in ", output_file)


if __name__=="__main__":
    main()