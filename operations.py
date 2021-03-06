import numpy as np
import math
from math import inf

from numpy.core.fromnumeric import sort

two_opt_swap = lambda r,i,k: np.concatenate((r[0:i+1],r[k:i:-1],r[k+1:len(r)]))


def get_abc(xs, ys):
    ''' Returns coeffs of the equations of each line:
        ax + by + c = 0. 
    '''
    a = ys[:-1] - ys[1:]
    b = xs[1:] - xs[:-1]
    c = xs[:-1]*ys[1:] - xs[1:]*ys[:-1]

    return a,b,c


def intersect_x(a1, a2, b1, b2, c1, c2):
    ''' Returns have these lines intersection, 
                x-coord of intersection point of two lines.
    '''
    if (b1*a2-b2*a1)!=0:
        y_int = (c2*a1-c1*a2)/(b1*a2-b2*a1)  
        x_int = (-c1-b1*y_int)/a1
        
        return True, x_int
    else:
        return False, 0


def ray_tracing(data, sorted_data, n):
    ''' Returns how many times ray originating from 
        the centroid and parallel to the x axis
        intersects contour made from data sorted by sorted_data.

        This function is using to understand if centroid is in contour.
        It means that contour is single otherwise not. 
    '''

    res = 0
    
    center = [np.sum(data[:, 0])/n,np.sum(data[:, 1])/n]
    max_x = np.max(data[:,0])

    n = len(sorted_data)
    xs = data[sorted_data,0]
    ys = data[sorted_data,1]
    xs = np.append(xs,xs[0])
    ys = np.append(ys,ys[0])
    a, b, c = get_abc(xs, ys)
    
    for i in range (n-1):
        if a[i]!=0:        
            x_int = (-c[i]-b[i]*center[1])/a[i]
            res += 1 if (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                        (center[0] <= x_int <= max_x)) else 0
        else:
            res += 1 if ys[i]==center[1] else 0

    return res


def count_loops(data, sorted_data, n, show=False):
    ''' Returns number of intersections in data sorted 
        in 'sorted_data' order if parameter show is False
        Otherwise returns arrays with indexes of sorted_data 
        where contour intersects.
    '''

    res = 0
    sorted_data = np.append(sorted_data,sorted_data[0])
    xs = data[sorted_data,0]
    ys = data[sorted_data,1]
    n = len(xs)
    beg = []
    end = []
    a,b,c = get_abc(xs, ys)
    for i in range (n-3):
        # consider non-adjacent edges
        for j in range (i+2, n-1):
            # don't consider last edge while considering first edge 
            # because they're adjacent
            if i==0 and j==n-2:
                continue
            # don't consider edge that goes after degenerated to a point edge
            if j==i+2 and (a[j-1]==0 and b[j-1]==0):
                continue
            not_null, x_int = intersect_x(a[i], a[j], b[i], b[j], c[i], c[j])

            if not_null and (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                                    min(xs[j], xs[j+1]) <= x_int <= max(xs[j], xs[j+1])):
                res+=1        
                beg.append(i)
                end.append(j)

    if show:
        return beg, end
    return res


def no_loops(data, sorted_data, n):
    ''' Returns data without loops using algorithm 
        based on 2-opt swapping mechanism. 
    '''

    x = data[:,0]
    y = data[:,1]
    new = sorted_data.copy()

    # Delete non-unique data from sorted array 'new', save it to 'non_uni'
    uni, inv, counts = np.unique(data[sorted_data], return_inverse=True, return_counts=True, axis=0)
    non_uni = []
    non_cons_points = []
    for i in range (len(uni)):
        if counts[i]>1:
            rep = np.where(inv==i)[0]
            non_uni.append(sorted_data[rep])
            non_cons_points.append(rep[1:])
    new = np.delete(new,non_cons_points)
    n = len(new) 
    new = np.append(new,new[0])

    # find intersections and delete them using 'two_opt_swap' function
    flag = True
    it = 0
    # normalize = lambda x: x / np.linalg.norm(x)
    # de = np.array([])
    while flag and it<50:
        it+=1
        flag = False
        
        xs = x[new]
        ys = y[new]
        
        a, b, c = get_abc(xs, ys)

        for i in range (n-2):
            # consider non-adjacent edges
            for j in range (i+2, n):
                # don't consider last edge while considering first edge 
                # because they're adjacent
                if i==0 and j==n-1:
                    continue
                
                # find intersection of two edges
                not_null, x_int = intersect_x(a[i], a[j], b[i], b[j], c[i], c[j])
                # check where the intersection lies
                if not_null and (min(xs[i], xs[i+1]) <= x_int <= max(xs[i], xs[i+1]) and
                                 min(xs[j], xs[j+1]) <= x_int <= max(xs[j], xs[j+1])):
                    flag = True
                    new = two_opt_swap(new,i,j)
                    break
            else:
                continue
            break

        # Not quite a successful attempt to remove deadlocks.
        # It removed dead ends, but because of this, loops appeared
        # And took more time to find contour and built contour less than 1% shorter 
        # for imporoved NN and N21 algorithms and less than 1% longer 
        # for improved CH algorithm.
        # Because of this this part was commented.

        # else:
        #     vec = np.vstack((-a, b)).T
        #     norm_vec = np.array(list(map(normalize, vec)))
        #     dotprod = norm_vec[:-1,0]*norm_vec[1:,0] + norm_vec[:-1,1]*norm_vec[1:,1]
        #     de1 = np.where(dotprod[1:]<-1+1e-3)[0]
        #     if len(de)==len(de1):
        #         if np.all(de==de1):
        #             break
        #     de = de1
        #     for ind in de:
        #         new = two_opt_swap(new,ind,ind+2)
        #         flag = True

    # add non-unique data after similar data
    for i in range (len(non_uni)):
        ind = np.where(new==non_uni[i][0])[0]
        new = np.insert(new,ind,non_uni[i][1:])

    # delete closing edge    
    new = np.delete(new,-1)
    return new


def sort_nn(data, sim, n):
    ''' Returns data sorted by the nearest neighbour algorithm. '''
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
 
    return np.array(sorted_data)     


def sort_nn_md(data, sim, n, return_md=False):
    ''' Returns data sorted by the nearest neighbour algorithm 
        which is single contour. 
    '''
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

    sorted_data = np.array(sorted_data)
    if return_md:
        missed_data = np.arange(n)
        md = np.delete(missed_data, missed_data[sorted_data])
        return sorted_data, md
    
    else:
        return sorted_data


def sort_nn_loops(data, sim, n):
    ''' Returns data sorted by the nearest neighbour algorithm
        and missed data appending algorithm. 
    '''    

    sorted_data, md = sort_nn_md(data, sim, n, return_md=True)
    sorted_data = append_md(sorted_data, md, sim)
    
    return sorted_data


def sort_nn_no_loops(data, sim, n):
    ''' Returns data sorted by the nearest neighbour algorithm,
        missed data appending algorithm and
        loops removal algorithm. 
    '''    

    sorted_data = sort_nn_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


def append_md(sorted_data, md, sim):
    ''' Returns sorted_data to which md was added 
        using nearest neighbour concept. 
    '''

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


def sort_21_loops(data, sim, n):
    ''' Returns data sorted by algorithm based on 
        inserting the second contour into the first one. 
    '''

    mid = np.argmin(sim[0][4:])
    sorted_data = np.arange(mid)
    md = np.delete(np.arange(n), sorted_data)

    sorted_data = append_md(sorted_data, md, sim)

    return sorted_data


def sort_21_no_loops(data, sim, n):
    ''' Returns data sorted by algorithm based on 
        inserting the second contour into the first one 
        and loop removal algorithm.
    '''

    sorted_data = sort_21_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


def clockwiseangle_and_distance(point, origin, refvec):
    ''' Returns polar coordinates of point with center = origin
        and axis = refvec. 
    '''

    vec = [point[0]-origin[0], point[1]-origin[1]]
    lenvec = math.hypot(vec[0], vec[1])
    
    if lenvec == 0:
        return -math.pi, 0
    
    normalized = [vec[0]/lenvec, vec[1]/lenvec]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    
    if angle < 0:
        return 2*math.pi+angle, lenvec
    
    return [angle, lenvec]


def sort_angle(data, sim, n):
    ''' Returns data sorted by polar coordinates. '''
    origin = [np.sum(data[:, 0])/n,np.sum(data[:, 1])/n]
    refvec = [0,1]
    # find polar coordinates
    angle = np.array([clockwiseangle_and_distance(point, origin, refvec) for point in data])
    # sort by polar coordinates
    sort = np.lexsort((angle[:,1],angle[:,0]))

    return sort


def sort_best(data, sim, n):
    ''' Returns data sorted by NN with no loops algorithm
        if it works correctly 
        otherwise data sorted by adding 2nd contour into 1st with no loops algorithm
        if it works correctly. 
    '''

    res = sort_nn_no_loops(data, sim, n)
    if count_loops(data, res, n)>0:
        res = sort_21_no_loops(data, sim, n)

    return res


def rotate(A,B,C):
    ''' Returns where is point C for vector AB 
        if returns negative number - C is right
        otherwise left.
    '''
    return (B[0]-A[0])*(C[1]-B[1])-(B[1]-A[1])*(C[0]-B[0])


def jarvismarch(data, sim, n):
    ''' Returns minimum convex polygon
        using Jarvis algorithm. 
    '''
    
    points = np.arange(n)
    ind_min_x = np.argmin(data[:,0])
    hull = [ind_min_x]
    points = np.delete(points,ind_min_x)
    points = np.append(points,ind_min_x)

    while True:
        right = 0
        for i in range(1,len(points)):
            if rotate(data[hull[-1]],data[points[right]],data[points[i]])<0:
                right = i
        if points[right]==hull[0]: 
            break
        else:
            hull = np.append(hull,points[right])
            points = np.delete(points,right)
    return hull


def distance_ch(p1,p2,p3):
    ''' Returns distance from point p3 to line p1,p2. '''
    d = np.cross(p2-p1,p3-p1)/np.linalg.norm(p2-p1)
    return d


def distance_ch_vert(d,a):
    ''' Returns hypotenuse of a right triangle with legs d and a. '''
    return np.sqrt(d**2+a**2)


def sort_ch_loops(data, sim, n):
    ''' Returns indexes of data sorted using algorithm 
        based on adding point into Convex Hull. 
    '''

    # find hull
    hull = jarvismarch(data, sim, n)
    nearest_line = [[i] for i in hull]
    distances = [[0] for i in hull]
    hull = np.append(hull, hull[0])
    
    # find nearest edge of convex hull and distance to it for other points
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
    
    # sort points by the distance of their projection from 
    # the first point of the nearest edge
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
    sorted_data = np.array(sorted_data)
    return sorted_data


def sort_ch_no_loops(data, sim, n):
    ''' Returns indexes of data sorted using algorithm 
        based on adding point into Convex Hull and 
        loop removal algorithm. 
    ''' 

    sorted_data = sort_ch_loops(data, sim, n)
    new = no_loops(data, sorted_data, n)

    return new


def two_opt(cities,sim,n,improvement_threshold=0.001): 
    ''' Returns indexes of data sorted using 2-opt algorithm. '''

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
    ''' Returns length of contour made of data sorted 
        by indexes in parameter sorted. 
    '''

    path = 0
    for i in range (n-1):
        path += sim[sorted[i]][sorted[i+1]]
    path += sim[sorted[-1]][sorted[0]]
  
    return path

