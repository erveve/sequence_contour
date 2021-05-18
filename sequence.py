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
             "ch": sort_ch,
             "2-opt": two_opt,
             "ui": sort_union_intersection}

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
        res = count_loops(self.data, self.contour, self.n)
        if return_num:
            return res>0, res
        else:
            return res>0

    
    def have_missed_data(self, return_num=False):
        # missed_data = np.arange(self.n)
        md = np.delete(np.arange(self.n), self.contour)
        if return_num:
            return len(md)>0, md
        else:
            return len(md)>0

    
    def is_data_unique(self, return_not_unique=False):
        _, indexes = np.unique(self.data, return_index=True, axis=0)
        res = np.delete(np.arange(self.n), indexes)
        # set(range(self.n))-set(indexes)
        if return_not_unique:
            return len(res)==0, res
        else:
            return len(res)==0


    def is_contour(self, return_problems=False):
        loop = self.have_loops()
        md = self.have_missed_data()
        unique_data = self.is_data_unique()
        res = True if not loop and not md and unique_data else False
        if return_problems:
            res_str = loop*'Have loops. '+md*'Have missed data. '+(not unique_data)*'Data is not unique.'
            return res, res_str
        else:
            return res


    def get_data_len(self):
        return self.n


    def get_contour_len(self):
        return calc_path(self.data, self.sim, self.contour, self.n)


    def sorted(self, key="nn_no_loops"):
        sort_ind = sort_dict[key](self.data, self.sim, self.n)
        sorted_data = Sequence(self.data, sort_ind)
        return sorted_data


    def show_contour(self, save_fig=False, fig_name='figure'):
        have_md, md = self.have_missed_data(return_num=True)

        plt.figure(num)
        if have_md:
            plt.scatter(self.data[md,0], self.data[md,1], c='g', s=30)
        plt.scatter(self.data[self.contour,0],self.data[self.contour,1], c='b', s=30)
        plt.plot(self.data[self.contour,0],self.data[self.contour,1],'r')
        if fig_name=='polar':
            plt.scatter(np.sum(self.data[:, 0])/self.n,np.sum(self.data[:, 1])/self.n, c='y', s=30)
        if save_fig:
            plt.savefig("images/"+fig_name+".png")
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
    