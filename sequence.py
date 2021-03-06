from matplotlib import pyplot as plt
from operations import *
import os

global sort_dict
sort_dict = {"nn": sort_nn,
             "nn_md": sort_nn_md,
             "nn_loops": sort_nn_loops,
             "nn_no_loops": sort_nn_no_loops,
             "nn_21_loops": sort_21_loops,
             "nn_21_no_loops": sort_21_no_loops,
             "polar": sort_angle,
             "ch": jarvismarch,
             "ch_loops": sort_ch_loops,
             "ch_no_loops": sort_ch_no_loops,
             "2-opt": two_opt,
             "best": sort_best}
    

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
        # print(res)
        if return_num:
            return res>0, res
        else:
            return res>0

    
    def have_missed_data(self, return_num=False):
        md = np.delete(np.arange(self.n), self.contour)
        if return_num:
            return len(md)>0, md
        else:
            return len(md)>0

    
    def is_data_unique(self, return_not_unique=False):
        _, indexes = np.unique(self.data, return_index=True, axis=0)
        res = np.delete(np.arange(self.n), indexes)
        if return_not_unique:
            return len(res)==0, res
        else:
            return len(res)==0

    
    def is_single_contour(self):
        intersections = ray_tracing(self.data, self.contour, self.n)
        return (intersections % 2)==1 

    
    def is_contour(self, return_all=False):
        loop = self.have_loops()
        md = self.have_missed_data()
        sc = self.is_single_contour()

        res = True if not loop and not md and sc else False

        if return_all:
            res_all = {'Missed data': md, 
                   'Loops': loop,
                   'Single contour': sc, 
                   'Contour': res}
            return res_all
        else:
            return res


    def get_data_len(self):
        return self.n


    def get_contour_len(self):
        return calc_path(self.data, self.sim, self.contour, self.n)


    def sorted(self, key="best"):
        sort_ind = sort_dict[key](self.data, self.sim, self.n)
        sorted_data = Sequence(self.data, sort_ind)
        return sorted_data


    def show_contour(self, save_fig=False, fig_name='figure'):
        have_md, md = self.have_missed_data(return_num=True)
        beg, end = count_loops(self.data, self.contour, self.n, show=True)

        plt.figure(num)
        if have_md:
            plt.scatter(self.data[md,0], self.data[md,1], c='g', s=30)
        plt.scatter(self.data[self.contour,0],self.data[self.contour,1], c='b', s=30)
        plt.plot(self.data[self.contour,0],self.data[self.contour,1],'r')
        for i in range (len(beg)):
            plt.plot(self.data[self.contour[beg[i]:beg[i]+2],0], 
                        self.data[self.contour[beg[i]:beg[i]+2],1],'y')
            plt.plot(self.data[self.contour[end[i]:end[i]+2],0], 
                        self.data[self.contour[end[i]:end[i]+2],1],'y')
        
        if fig_name=='polar':
            plt.scatter(np.sum(self.data[:, 0])/self.n,np.sum(self.data[:, 1])/self.n, c='y', s=30)
        if save_fig:
            plt.savefig("images/"+fig_name+".png")
        else:
            plt.show()        


if __name__=="__main__":
    # save pictures for diploma
    datas = []
    for subdir, dirs, files in os.walk('datasets'):
        for file in files:
            filepath = subdir + os.sep + file
            if filepath.endswith(".dat"):
                datas.append(np.loadtxt(filepath))
    n = len(datas)
    
    seq = Sequence(datas[0]) 
    global num
    num = 0
    seq.show_contour(save_fig=True)
    for method in sort_dict:
        num+=1
        sort_seq = seq.sorted(key=method)
        sort_seq.show_contour(save_fig=True, fig_name=method)
        # sort_seq.show_contour()
    