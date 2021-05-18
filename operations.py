import numpy as np
import math
from math import inf

from shapely.geometry import Point, Polygon

two_opt_swap = lambda r,i,k: np.concatenate((r[0:i+1],r[k:i:-1],r[k+1:len(r)]))

def count_loops(data, sorted_data, n):
    res = 0

    xs = data[sorted_data,0]
    ys = data[sorted_data,1]

    a = [ys[i]-ys[i+1] for i in range (n-1)]
    b = [xs[i+1]-xs[i] for i in range (n-1)]
    c = [xs[i]*ys[i+1] - xs[i+1]*ys[i] for i in range (n-1)]
    for i in range (n-3):
        for j in range (i+2, n-1):
            if (b[i]*a[j]-b[j]*a[i]!=0):
                y_int = (c[j]*a[i]-c[i]*a[j])/(b[i]*a[j]-b[j]*a[i])
                x_int = (-c[i]-b[i]*y_int)/a[i]
                res += 1 if (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                             min(xs[j], xs[j+1]) <= x_int <= max(xs[j], xs[j+1])) else 0
                    
    return res


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
            dist = sim[i][sorted_data[j]] + sim[i][sorted_data[j+1]]
            # ij = np.linalg.norm(data[i]-data[sorted_data[j]])
            # ij1 = np.linalg.norm(data[i]-data[sorted_data[j+1]])
            if dist<tmp:
                tmp = dist
                ind = j 
        if sim[i][sorted_data[-1]] + sim[i][sorted_data[0]] < tmp:
        # if np.linalg.norm(data[i]-data[sorted_data[j]])+np.linalg.norm(data[i]-data[0])<tmp:
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


def rotate(A,B,C):
    ''' returns where is point C for vector AB 
        if returns negative number - C is right
        otherwise left'''
    return (B[0]-A[0])*(C[1]-B[1])-(B[1]-A[1])*(C[0]-B[0])


def grahamscan(A, n):
    ''' returns minimum convex polygon
        using Graham's algorithm '''
    P = [i for i in range(n)] # list with points' numbers
    # find point with min x
    for i in range(1,n):
        if A[P[i]][0]<A[P[0]][0]: # if point P[i] is left from point P[0]
            P[i], P[0] = P[0], P[i] # swap their numbers
    
    # sort points против часовой стрелки
    for i in range(2,n):
        j = i
        while j>1 and (rotate(A[P[0]],A[P[j-1]],A[P[j]])<0): 
            P[j], P[j-1] = P[j-1], P[j]
            j -= 1

    # cut edges
    S = [P[0],P[1]] # first points are in min 
    for i in range(2,n):
        while rotate(A[S[-2]],A[S[-1]],A[P[i]])<0:
            del S[-1] # pop(S)
        S.append(P[i]) # push(S,P[i])
    
    return S


def jarvismarch(A, n):
    ''' returns minimum convex polygon
        using Jarvis algorithm '''
    P = [i for i in range(n)] # list with points' numbers
    for i in range(1,n):
        if A[P[i]][0]<A[P[0]][0]: 
            P[i], P[0] = P[0], P[i] 
    
    H = [P[0]]
    del P[0]
    P.append(H[0])

    while True:
        right = 0
        for i in range(1,len(P)):
            if rotate(A[H[-1]],A[P[right]],A[P[i]])<0:
                right = i
        if P[right]==H[0]: 
            break
        else:
            H.append(P[right])
            del P[right]
    return H


def distance_ch(p1,p2,p3):
    d = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
    return d


def distance_ch_vert(d,a):
    return np.sqrt(d**2+a**2)


def sort_ch(data, sim, n):
    hull = jarvismarch(data, n)
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
    

def calc_path(data, sim, sorted, n):
    path = 0
    for i in range (n-1):
        path += sim[sorted[i]][sorted[i+1]]
        # path += np.linalg.norm(data[sorted[i+1]]-data[sorted[i]])
    path += sim[sorted[-1]][sorted[0]]
    # path += np.linalg.norm(data[sorted[-1]]-data[sorted[0]])
  
    return path


def two_opt(cities,sim,n,improvement_threshold=0.001): 
    route = np.arange(n) 
    improvement_factor = 1 
    # the distance of the initial path.
    best_distance = calc_path(cities,sim,route,n)
    while improvement_factor > improvement_threshold: 
        distance_to_beat = best_distance 
        for swap_first in range(1,n-2): 
            for swap_last in range(swap_first+1,n): 
                new_route = two_opt_swap(route,swap_first,swap_last) 
                new_distance = calc_path(cities,sim,new_route,n) 
                if new_distance < best_distance: 
                    route = new_route 
                    best_distance = new_distance 
        improvement_factor = 1 - best_distance/distance_to_beat 
    return route 


def split_2_contours(data,sim,n):
    # split data in two contours
    mid = np.argmin(sim[0][4:])-10
    mid = 5
    while mid<n-10:
        data1 = data[:mid]
        data2 = data[mid:]

        # arr1 = [Point(data1[i][0],data1[i][1]) for i in range (mid)]
        arr1 = [Point(data1[i][0],data1[i][1]) for i in range (len(data1))]
        arr1.append(arr1[0])
        # arr2 = [Point(data2[i][0],data2[i][1]) for i in range (mid)]
        arr2 = [Point(data2[i][0],data2[i][1]) for i in range (len(data2))]
        arr2.append(arr2[0])
        lr1 = Polygon(arr1)
        lr2 = Polygon(arr2)

        if lr1.is_valid and lr2.is_valid:
            break
        mid += 1

    if mid==n-10:
        mid = np.argmin(sim[0][4:])+2
        data1 = data[:mid]
        data2 = data[mid:]

        # arr1 = [Point(data1[i][0],data1[i][1]) for i in range (mid)]
        arr1 = [Point(data1[i][0],data1[i][1]) for i in range (len(data1))]
        arr1.append(arr1[0])
        # arr2 = [Point(data2[i][0],data2[i][1]) for i in range (n-mid)]
        arr2 = [Point(data2[i][0],data2[i][1]) for i in range (len(data2))]
        arr2.append(arr2[0])
        lr1 = Polygon(arr1)
        lr2 = Polygon(arr2)
    
    if not lr1.is_valid:
        mid = np.argmin(sim[0][4:])+2
        data1 = data[no_loops(data[:mid],list(range(mid)))]
        # arr1 = [Point(data1[i][0],data1[i][1]) for i in range (mid)]
        arr1 = [Point(data1[i][0],data1[i][1]) for i in range (len(data1))]
        arr1.append(arr1[0])
        lr1 = Polygon(arr1)

    if not lr2.is_valid:
        mid = np.argmin(sim[0][4:])+2
        ind = no_loops(data[mid:],list(range(n-mid)))
        ind = np.array(ind)+mid
        data2 = data[ind]
        # arr2 = [Point(data2[i][0],data2[i][1]) for i in range (n-mid)]
        arr2 = [Point(data2[i][0],data2[i][1]) for i in range (len(data2))]
        arr2.append(arr2[0])
        lr2 = Polygon(arr2) 

    
    return lr1, lr2


def get_union_intersction(data,sim,n):
    lr1,lr2 = split_2_contours(data,sim,n)
    ux,uy = lr1.union(lr2).boundary.coords.xy
    union_data = np.array([-1]*len(ux))
    ix,iy = lr1.intersection(lr2).boundary.coords.xy
    intersect_data = np.array([-1]*len(ux))
    eps = 0.0001
    for i in range(n):
        for j in range(len(ux)):
            if abs(ux[j]-data[i,0])<=eps and abs(uy[j]-data[i,1])<=eps:
                union_data[j] = i
    #             break
        for j in range(len(ix)):
            if abs(ix[j]-data[i,0])<=eps and abs(iy[j]-data[i,1])<=eps:
                intersect_data[j] = i
    #             break

    union_data = union_data[union_data>-1]
    union_data = union_data[:-1]

    intersect_data = intersect_data[intersect_data>-1]
    intersect_data = intersect_data[:-1]
    return union_data, intersect_data

    
def sort_union_intersection(data,sim,n):
    union_data, intersect_data = get_union_intersction(data, sim, n)
    sort_ui = union_data
    intersect_data = np.append(intersect_data,list(set(range(n))-(set(union_data)|set(intersect_data))))
    intersect_data = intersect_data.astype(int)

    for el in intersect_data:
        ij = np.linalg.norm(data[el]-data[sort_ui[:-1]],axis=1)
        ij1 = np.linalg.norm(data[el]-data[sort_ui[1:]],axis=1)
        ij = ij + ij1
    #     ind = np.where(ij==np.amin(ij))[0][0]
        ind = np.argmin(ij)
        sort_ui = np.insert(sort_ui,ind+1,el)
        
    return sort_ui    