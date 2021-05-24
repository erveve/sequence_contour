import numpy as np
import math
from math import inf

two_opt_swap = lambda r,i,k: np.concatenate((r[0:i+1],r[k:i:-1],r[k+1:len(r)]))


def get_abc(xs, ys, n):
    ''' Returns coeffs of the equations of each line:
        ax + by + c = 0 '''
    a = ys[:-1] - ys[1:]
    b = xs[1:] - xs[:-1]
    c = xs[:-1]*ys[1:] - xs[1:]*ys[:-1]

    return a,b,c


def intersect_x(a1, a2, b1, b2, c1, c2):
    ''' Returns have these lines intersection, 
                x-coord of intersection point of two lines'''
    if (b1*a2-b2*a1)!=0:
        y_int = (c2*a1-c1*a2)/(b1*a2-b2*a1)  
        x_int = (-c1-b1*y_int)/a1
        
        return True, x_int
    else:
        return False, 0


def ray_tracing(data, sorted_data, n):
    res = 0
    
    center = [np.sum(data[:, 0])/n,np.sum(data[:, 1])/n]
    max_x = np.max(data[:,0])

    n = len(sorted_data)
    xs = data[sorted_data,0]
    ys = data[sorted_data,1]
    xs = np.append(xs,xs[0])
    ys = np.append(ys,ys[0])
    a, b, c = get_abc(xs, ys, n)
    
    for i in range (n-1):
        if a[i]!=0:        
            x_int = (-c[i]-b[i]*center[1])/a[i]
            res += 1 if (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                        (center[0] <= x_int <= max_x)) else 0
        else:
            res += 1 if ys[i]==center[1] else 0

    return res


# def get_dead_ends(data, sorted_data):
#     _, indices = np.unique(data[sorted_data], return_inverse=True, axis=0)
#     res = np.append(indices[:2]==indices[-2:], indices[2:]==indices[:-2])
#     return np.where(res)[0]


def count_loops(data, sorted_data, n):
    res = 0
    de = 0

    xs = data[sorted_data,0]
    ys = data[sorted_data,1]
    n = len(xs)

    a,b,c = get_abc(xs, ys, n)
    for i in range (n-3):
        for j in range (i+2, n-1):
            not_null, x_int = intersect_x(a[i], a[j], b[i], b[j], c[i], c[j])
            res += 1 if not_null and (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                                      min(xs[j], xs[j+1]) <= x_int <= max(xs[j], xs[j+1])) else 0

            de += 1 if not not_null and c[i]*a[j]==c[j]*a[i] and not(xs[i]==xs[j] and ys[i]==ys[j]) else 0

    return res, de


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

        a, b, c = get_abc(xs, ys, n)

        for i in range (n-3):
            for j in range (i+2, n-1):
                not_null, x_int = intersect_x(a[i], a[j], b[i], b[j], c[i], c[j])
                if not_null and (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
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

            if dist<tmp:
                tmp = dist
                ind = j 
        if sim[i][sorted_data[-1]] + sim[i][sorted_data[0]] < tmp:
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


def sort_best(data, sim, n):
    res = [sort_nn_no_loops(data, sim, n)]
    cnt, _ = count_loops(data, res[0], n)
    cnt_loops = [cnt]
    if cnt_loops[-1]>0:
        res.append(sort_21_no_loops(data, sim, n))
        cnt, _ = count_loops(data, res[-1], n)
        cnt_loops.append(cnt)
        if cnt_loops[-1]>0:
            res.append(sort_ch_no_loops(data, sim, n))
            cnt, _ = count_loops(data, res[-1], n)
            cnt_loops.append(cnt)

    return res[np.argmin(np.array(cnt_loops))]


def rotate(A,B,C):
    ''' returns where is point C for vector AB 
        if returns negative number - C is right
        otherwise left'''
    return (B[0]-A[0])*(C[1]-B[1])-(B[1]-A[1])*(C[0]-B[0])


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


def sort_ch_loops(data, sim, n):
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


def sort_ch_no_loops(data, sim, n):
    sorted_data = sort_ch_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


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


def calc_path(data, sim, sorted, n):
    path = 0
    for i in range (n-1):
        path += sim[sorted[i]][sorted[i+1]]
    path += sim[sorted[-1]][sorted[0]]
  
    return path

