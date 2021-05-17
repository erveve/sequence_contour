import numpy as np
import math
from math import inf
from scipy.spatial import ConvexHull

two_opt_swap = lambda r,i,k: np.concatenate((r[0:i+1],r[k:i:-1],r[k+1:len(r)]))

def no_loops(data, sorted_data, n):
    x = data[:,0]
    y = data[:,1]
    new = sorted_data.copy()
    flag = True
    it = 0
    while flag and it<50:
        it+=1
        flag = False
        
        xs = x[new]
        ys = y[new]

        a = [ys[i]-ys[i+1] for i in range (n-1)]
        b = [xs[i+1]-xs[i] for i in range (n-1)]
        c = [xs[i]*ys[i+1] - xs[i+1]*ys[i] for i in range (n-1)]
        for i in range (n-3):
            for j in range (i+2, n-1):
                if (b[i]*a[j]-b[j]*a[i]!=0):
                    y_int = (c[j]*a[i]-c[i]*a[j])/(b[i]*a[j]-b[j]*a[i])
                    x_int = (-c[i]-b[i]*y_int)/a[i]
                    if (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                        min(xs[j], xs[j+1]) <= x_int <= max(xs[j], xs[j+1])):
                        flag = True
                        new = two_opt_swap(new,i,j)
                        break
            else:
                continue
            break

    return new


def sort_nn(data, sim, n):
    sorted_data = [0]
    while True:
        tmp = inf
        ind = -1    
        for i in range (n):
            if i not in sorted_data:
                if sim[sorted_data[-1],i]<tmp:
                    tmp = sim[sorted_data[-1],i]
                    ind = i
        if ind>0:
            sorted_data.append(ind)
        else:
            break        
    return sorted_data


def sort_nn_md(data, sim, n, return_md=False):
    sorted_data = [0]
    while True:
        tmp = inf
        ind = -1    
        for i in range (n):
            if i not in sorted_data:
                if sim[sorted_data[-1],i]<tmp:
                    tmp = sim[sorted_data[-1],i]
                    ind = i
        if len(sorted_data)>10 and sim[sorted_data[-1],0]<tmp:
            break
        if ind>0:
            sorted_data.append(ind)
        else:
            break    
            
    if return_md:
        missed_data = np.arange(n)
        md = np.delete(missed_data, missed_data[sorted_data])
        return sorted_data, md
    
    else:
        return sorted_data


def sort_nn_loops(data, sim, n):
    sorted_data, md = sort_nn_md(data, sim, n, return_md=True)

    for i in md:
        tmp = inf
        ind = -1
        for j in range (len(sorted_data)-1):
            ij = np.linalg.norm(data[i]-data[sorted_data[j]])
            ij1 = np.linalg.norm(data[i]-data[sorted_data[j+1]])
            if (ij + ij1)<tmp:
                tmp = ij +ij1
                ind = j                                                                         
        if np.linalg.norm(data[i]-data[sorted_data[j]])+np.linalg.norm(data[i]-data[0])<tmp:
            sorted_data.append(i)
        else:
            sorted_data.insert(ind+1,i)
    
    return sorted_data


def sort_nn_no_loops(data, sim, n):
    sorted_data = sort_nn_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


def sort_21_loops(data, sim, n):
    mid = np.argmin(sim[0][4:])
    data1 = data[:mid]
    data2 = data[mid:]

    sorted_data = np.arange(mid)
    missed_data = np.arange(n)
    md = np.delete(missed_data, missed_data[sorted_data])

    for i in md:
        tmp = inf
        ind = -1
        for j in range (len(sorted_data)-1):
            ij = np.linalg.norm(data[i]-data[sorted_data[j]])
            ij1 = np.linalg.norm(data[i]-data[sorted_data[j+1]])
            if (ij + ij1)<tmp:
                tmp = ij +ij1
                ind = j                                                                         
        if np.linalg.norm(data[i]-data[sorted_data[j]])+np.linalg.norm(data[i]-data[0])<tmp:
            sorted_data = np.append(sorted_data,i)
        else:
            sorted_data = np.insert(sorted_data,ind+1,i)

    return sorted_data


def sort_21_no_loops(data, sim, n):
    sorted_data = sort_21_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


def clockwiseangle_and_distance(point, origin, refvec):
    vector = [point[0]-origin[0], point[1]-origin[1]]
    lenvector = math.hypot(vector[0], vector[1])
    
    if lenvector == 0:
        return -math.pi, 0
    
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    
    if angle < 0:
        return 2*math.pi+angle, lenvector
    
    return angle, lenvector


def sort_angle(data, sim, n):
    origin = [np.sum(data[:, 0])/n,np.sum(data[:, 1])/n]
    refvec = [0,1]
    angle = np.array([clockwiseangle_and_distance(point, origin, refvec) for point in data])
    sort = np.argsort(angle[:,0])

    return sort


# def convex_hull(data, n):


def distance_ch(p1,p2,p3):
    d = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
    return d


def distance_ch_vert(d,a):
    return np.sqrt(d**2+a**2)


def sort_ch(data, sim, n):
    hull = ConvexHull(data)
    hull = hull.vertices
    nearest_line = [[i] for i in hull]
    distances = [[0] for i in hull]
    hull = np.append(hull, hull[0])
    
    for i, point in enumerate(data):
        if i in hull:
            continue
        else:
            a = inf
            ind = 0
            for j in range (len(hull)-1):
                d = distance_ch(data[hull[j]],data[hull[j+1]],point)
                if abs(d)<abs(a):
                    a = d
                    ind = j
            nearest_line[ind].append(i)
            distances[ind].append(a)
    
    sorted_data = []
    for i in range(len(nearest_line)):
        c = nearest_line[i]
        a = distances[i]
        b =[[i] for i in range(len(c))]
        b[0].append(0)
        for j in range(1, len(c)):
            b[j].append(distance_ch_vert(sim[c[0],c[j]],a[j]))
        b = sorted(b,key=lambda dist:dist[1])
        
        for j in range(len(c)):
            sorted_data.append(c[b[j][0]])
    
    return sorted_data
    

def calc_path(data, sorted, n):
    path = 0
    for i in range (n-1):
        path += np.linalg.norm(data[sorted[i+1]]-data[sorted[i]])
    path += np.linalg.norm(data[sorted[-1]]-data[sorted[0]])
  
    return path