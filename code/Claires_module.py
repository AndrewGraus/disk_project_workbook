#Here I'm going to put Claire's module

pnts = 100

DES_cone_len = 1.0 #Satellite[lstn][5]/1000.
SDSS_cone_len = 1.0 #Satellite[lstn][4]/1000.  #t2 3

DES_cone_size = 0.7259
SDSS_coneS_size = 0.7259  #t2
SDSS_coneB_size = 0.9929   #t3

#corr_factors = np.zeros((len(data_files),pnts))
#count_mas = np.zeros((len(data_files),pnts))

radius_bins = np.linspace(0.0,0.1,101)

h = 0.6751

def cart2sph(x,y,z):
    XsqPlusYsq = np.square(x) + np.square(y)
    r = np.sqrt(XsqPlusYsq + np.square(z))         # r
    elev = np.arctan2(z,np.sqrt(XsqPlusYsq))       # phi
    az = np.arctan2(y,x)                           # theta
    return r, elev, az

def calc_open_angle(phi_x,the_x,phi_t,the_t):
    return np.arccos(np.cos(phi_t) * np.cos(phi_x) + np.sin(phi_t) * np.sin(phi_x) * np.cos(the_x - the_t))

def calc_corr_factors(Data_X, Data_Y, Data_Z, Data_M, Data_Rvir, pnts, MW, Andr):
    # Rvir converted in Mpc (was kpc) 
    # In order to ignore this I think I just need to set this to zero
    #RO2 = Data_Rvir[Andr]/1000.
    RO2 = 0.0
    

    # in Mpc
    z0_cens = np.vstack([Data_X,Data_Y,Data_Z]).T
    pair_dist = z0_cens[Andr] - z0_cens[MW]   # First two are the hosts (M31 & MW)
    coord_dist = z0_cens[2:] - z0_cens[MW]
    
    dis_g = np.linalg.norm(pair_dist)
    distances_from_center = np.linalg.norm(coord_dist,axis=1)[newaxis]
    
    angle_BG = np.arctan(RO2/dis_g)
    
    _, phi_x, the_x = cart2sph(coord_dist[:,0],coord_dist[:,1],coord_dist[:,2])
    _, phi_t, the_t = cart2sph(-pair_dist[0],-pair_dist[1],-pair_dist[2])

    ang_dis = calc_open_angle(phi_x,the_x,phi_t,the_t)

    Andr_gal = ang_dis <= angle_BG # Mask things that are in the same area as Andr.
    Andr_gal = Andr_gal[np.newaxis].T & np.ones(pnts,dtype = 'bool') # Make into N x pnts matrix

    #random
    #np.random.seed(777)
    u = np.random.uniform(0,1,pnts)
    v = np.random.uniform(0,1,pnts)

    the_t = 2. * np.pi * u
    phi_t = np.arccos(2. * v - 1)

    the_t2 = the_t - 0.09
    phi_t2 = phi_t + 1.18
    the_t3 = the_t + 3.1
    phi_t3 = phi_t + 1.72

    the_t2[ the_t2 < 0 ] += 2. * np.pi
    phi_t2[ phi_t2 > np.pi ] -= np.pi
    the_t3[ the_t3 > 2*np.pi ] -= 2*np.pi
    phi_t3[ phi_t3 > np.pi ] -= np.pi

    s = cart2sph(coord_dist[:,0],coord_dist[:,1],coord_dist[:,2])

    phi_x = (s[1])[np.newaxis].T
    the_x = (s[2])[np.newaxis].T
    r = s[0]

    dis = calc_open_angle(phi_x,the_x,phi_t,the_t)
    dis2 = calc_open_angle(phi_x,the_x,phi_t2,the_t2)
    dis3 = calc_open_angle(phi_x,the_x,phi_t3,the_t3)

    r_msk = (r[np.newaxis].T < np.ones(pnts)) & (~Andr_gal)

    r_SDSS = (r[np.newaxis].T < np.ones(pnts) * SDSS_cone_len) & r_msk
    s_SDSS = (dis2 <= SDSS_coneS_size) & r_SDSS
    b_SDSS = (dis3 <= SDSS_coneB_size) & r_SDSS

    in_SDSS = s_SDSS | b_SDSS
    in_DES = (r[np.newaxis].T < np.ones(pnts) * DES_cone_len) & (dis <= DES_cone_size) & r_msk

    in_cones = in_SDSS | in_DES
    
    #okay I believe if I just apply the below line mask to a matrix that is just the 
    #radii of each subhalo array stacked in a matrix 100 times.
    
    #for whatever reason the masked array doesn't really seem to work, and I really can't figure out why
    #so what I'm going to do is loop through each mask (weak) and then just go ahead and bin it up in a 
    #histogram since otherwise the different masks will lead to 
    
    stacked_distance_array_fast = np.broadcast_to(distances_from_center.T,(len(distances_from_center[0]),pnts))
    
    #note: the below line "works" but it gives a 1d array of all the values that pass ANY of the masks (not useful)
    #all_dist_in_cones = stacked_distance_array_fast[in_cones]
    
    #This doesn't work for reasons that escape me, it doesn't mask ANY values
    #distances_within_cones = np.ma.MaskedArray(stacked_distance_array_fast, mask = in_cones)
    
    #The below works but it uses a loop
    
    histogram_matrix = np.zeros((pnts,len(radius_bins)-1))
    
    for ii in range(len(stacked_distance_array_fast[0])):
        #pull out each individual mask
        indiv_mask = in_cones[:,ii]
        #pull out each of the distance arrays (note: I could just use the distances as in always using the
        #distances_from_center array, but I matricized it in case I ever figure out how to do this without
        #a loop)
        indiv_dist = stacked_distance_array_fast[:,ii]
        
        halos_within_cone = indiv_dist[indiv_mask]
        
        #now lets bin it with a histogram
        
        hist, bins = np.histogram(halos_within_cone,bins=radius_bins)
        histogram_matrix[ii] = np.cumsum(hist)
        
    count_mas = np.sum(in_cones[Data_M[2:] >= M_t], axis = 0 ) * 1.
    count_tot = np.sum(r_msk[Data_M[2:] >= M_t], axis = 0 ) * 1.
    
    return histogram_matrix
