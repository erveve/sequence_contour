import sys
import gmsh
from sequence import *

def load_data(input_file):
    ''' Load data from input_file '''
    data = np.loadtxt(input_file)
    return data


def sort_data(data):
    ''' Returns sorted data using best algorithm '''
    seq = Sequence(data) 
    sort_seq = seq.sorted(key='best')
    return sort_seq


def make_cad_file(sort_seq, output_file):
    ''' Makes CAD file from sort_seq and saves it in output_file '''
    model = gmsh.model
    factory = model.geo

    gmsh.initialize()

    model.add("contour")
    lc = 1e-5

    for i in range(sort_seq.n):
        factory.addPoint(sort_seq.data[i,0], sort_seq.data[i,1], 0, lc, i)
    
    for i in range(sort_seq.n-1):
        factory.addLine(sort_seq.contour[i], sort_seq.contour[i+1], i)
    factory.addLine(sort_seq.contour[-1], sort_seq.contour[0], sort_seq.n-1)

    factory.synchronize()

    gmsh.write(output_file)
    gmsh.finalize()


def main():

    if len(sys.argv)!=3:
        print("ERROR: Invalid number of arguments for main.py\nWrite: python main.py <input_file.dat> <output_file.geo_unrolled>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

        data = load_data(input_file)
        sort_seq = sort_data(data)
        make_cad_file(sort_seq, output_file)
        # np.savetxt(output_file, sort_seq.data[sort_seq.contour])


if __name__=="__main__":
    main()