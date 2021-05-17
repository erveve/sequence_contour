from ast import Num
import pandas as pd
from matplotlib import pyplot as plt
from operations import *

global sort_dict
sort_dict = {"nn": sort_nn,
             "nn_md": sort_nn_md,
             "nn_loops": sort_nn_loops,
             "nn_no_loops": sort_nn_no_loops,
             "nn_21_loops": sort_21_loops,
             "nn_21_no_loops": sort_21_no_loops,
             "polar": sort_angle,
             "ch": sort_ch}

class Sequence():
    def __init__(self, data, contour=[]):
        self.data = np.array(data)
        self.n = len(self.data)
        self.sim = np.array([[np.linalg.norm(self.data[i]-self.data[j]) for i in range (self.n)] for j in range (self.n)])
        if len(contour)<2:
            self.contour = np.array(range(self.n))
        else:
            self.contour = np.array(contour)

    
    def have_loops(self, return_num=False):
        return False

    
    def have_missed_data(self, return_num=False):
        if return_num:
            return False, 0
        return False

    
    def is_data_unique(self, return_not_unique=False):
        return True


    def is_contour(self):
        if not self.have_loops() and not self.have_missed_data:
            return True
        else:
            return False


    def get_contour_len(self):
        return len(self.contour)

    
    def sorted(self, key="nn_no_loops"):
        sort_ind = sort_dict[key](self.data, self.sim, self.n)
        sorted_data = Sequence(self.data, sort_ind)
        return sorted_data


    def show_contour(self, save_fig=False, fig_name='figure'):
        have_md, md = self.have_missed_data(return_num=True)
        if have_md:
            plt.scatter(self.data[md,0], self.data[md,1], c='g', s=30)
        plt.figure(num)
        plt.scatter(self.data[self.contour,0],self.data[self.contour,1], c='b', s=30)
        plt.plot(self.data[self.contour,0],self.data[self.contour,1],'r')
        if save_fig:
            plt.savefig(fig_name+".png")
        else:
            plt.show()        


if __name__=="__main__":
    data = pd.read_table("../data1/1.1.dat", sep=' ', header=None)
    data = data.to_numpy()
    seq = Sequence(data) 
    global num
    num = 0
    seq.show_contour(save_fig=True)
    for method in sort_dict:
        num+=1
        sort_seq = seq.sorted(key=method)
        sort_seq.show_contour(save_fig=True, fig_name=method)
    