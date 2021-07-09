import h5py, os
import numpy as np

#I'm doing something else here and grabbing the 
#from the merger tree files, there are less halos
#in the merger tree for a reason that I don't understand
#I'm hoping its because there very small halos are cut out

a_target = 1.0/(1.0+3.0) #want to find halos near z=3
h=0.6751

sim_type = 'dmo'

if sim_type == 'dmo':
    halo_dir = './dmo_files_full'
    snap_loc = 'Snapshot00152'
    output_dir = 'top_ten_files_dmo'
    tree_dir = 'dmo_files_full'
elif sim_type == 'disk':
    halo_dir = './disk_file_full'
    snap_loc = 'Snapshot000152'
    output_dir = 'top_ten_files'
    tree_dir = 'disk_files_full'
    rvir_param = 'rvir'

for data_file in os.listdir(halo_dir):
    if sim_type == 'dmo':
        halo_number = data_file.split('_')[2]
    else:
        halo_number = data_file.split('_')[1]

    print(data_file)
    
    #loading data
    f = h5py.File('./'+tree_dir+'/'+data_file)
    #print f.keys()

    snap = f[snap_loc]['HaloCatalog_RockstarMergerTree']
    Mvir = snap['Mvir'][:]
    Rvir = snap['Rvir'][:]
    cens = snap['Center'][:]
    Vmax = snap['Vmax'][:]
    #Rmax = snap['Rmax'][:]
    Vpeak = snap['Vpeak_TotalMass'][:]
    ids_from_halo_cat = snap['Orig_halo_id'][:]
    Tree_id = snap['TreeNumber'][:]
    Tree_root_id = snap['Tree_root_id'][:]
    scale = snap['Scale'][:]
    Vel = snap['Velocity'][:]

    #now I need the most massive halo at z = 0
    host_index = np.argmax(Mvir) #I'm assuming I can find the host by looking for the most massive thing in the entire tree (which should be the host at z=0
    host_tree_id = Tree_id[host_index] #This is the id of the host tree
    print('host_id: {}'.format(host_tree_id))
    #now grab the host
    host_tree = f['RockstarMergerTrees']['Tree_'+str(host_tree_id)]
    scale = host_tree['scale'][:]
    mmp = host_tree['mmp?'][:]
    Mvir_tree = host_tree['Mvir_all'][:]/h #Msun/h coverted to Msun
    x,y,z = host_tree['x'][:]*1000.0/h,host_tree['y'][:]*1000.0/h,host_tree['z'][:]*1000.0/h #cMpc/h converted to ckpc
    vx,vy,vz = host_tree['vx'][:],host_tree['vy'][:],host_tree['vz'][:]
    Tree_root_id = host_tree['Tree_root_ID'][:]
    #annoyingly sometimes rvir is Rvir and other times it's rvir so I guess I'll do a try
    try:
        rvir_tree = host_tree['rvir'][:]/h #ckpc/h converted to ckpc
    except KeyError:
        rvir_tree = host_tree['Rvir'][:]/h #ckpc/h converted to ckpc
    rs_tree = host_tree['rs'][:]/h #ckpc/h converted to ckpc

    #Now find the scale closest to 3
    unique_scale = np.unique(scale)
    scale_select = unique_scale[np.argmin(np.abs(unique_scale-a_target))]

    #everything from z=3 in the main tree
    #I don't think I need anything from the tree_root_id because they are all by defintion destroyed by z=0 
    #or they wouldn't be in the main branch
    print(scale_select)
    destroyed_three_mask = (scale==scale_select)
    
    print('halos at z=3 that are destroyed by z=0: {}'.format(np.sum(destroyed_three_mask)))

    destroyed_three_mass = Mvir_tree[destroyed_three_mask]
    destroyed_three_x, destroyed_three_y, destroyed_three_z = x[destroyed_three_mask]*scale_select,y[destroyed_three_mask]*scale_select,z[destroyed_three_mask]*scale_select
    destroyed_three_vx, destroyed_three_vy, destroyed_three_vz = vx[destroyed_three_mask],vy[destroyed_three_mask],vz[destroyed_three_mask]
    three_rvir = rvir_tree[destroyed_three_mask]*scale_select
    three_rs = rs_tree[destroyed_three_mask]*scale_select
    
    sorted_masses_three = np.sort(destroyed_three_mass)/1.0e11

    print('Most massive halo {}, second most massive {}'.format(sorted_masses_three[-1],sorted_masses_three[-2]))

    #Find the most massive halo
    host_halo = np.argmax(destroyed_three_mass)
    mmp_mass = destroyed_three_mass[host_halo]
    mmp_x,mmp_y,mmp_z = destroyed_three_x[host_halo],destroyed_three_y[host_halo],destroyed_three_z[host_halo]
    mmp_vx,mmp_vy,mmp_vz = destroyed_three_vx[host_halo],destroyed_three_vy[host_halo],destroyed_three_vz[host_halo]
    mmp_rvir = three_rvir[host_halo]
    mmp_rs = three_rs[host_halo]
    

    #now make sure it's a satellite at z=3
    x_diff,y_diff,z_diff = destroyed_three_x-mmp_x,destroyed_three_y-mmp_y,destroyed_three_z-mmp_z
    vx_diff,vy_diff,vz_diff = destroyed_three_vx-mmp_vx,destroyed_three_vy-mmp_vy,destroyed_three_vz-mmp_vz

    dist = np.sqrt(x_diff**2.0+y_diff**2.0+z_diff**2.0)
    dist_mask = (dist<mmp_rvir)&(dist>0.0)
    print('host Rvir at z=3: {} kpc'.format(mmp_rvir))
    print('destroyed satellites present at z=3: {}'.format(np.sum(dist_mask)))
    #everything that is a satellite
    destroyed_sat_three_mass = destroyed_three_mass[dist_mask]
    destroyed_sat_three_x, destroyed_sat_three_y, destroyed_sat_three_z = destroyed_three_x[dist_mask], destroyed_three_y[dist_mask], destroyed_three_z[dist_mask]
    destroyed_sat_three_vx, destroyed_sat_three_vy, destroyed_sat_three_vz = destroyed_three_vx[dist_mask], destroyed_three_vy[dist_mask], destroyed_three_vz[dist_mask]
    destroyed_sat_dist = dist[dist_mask]
    destroyed_sat_rvir = three_rvir[dist_mask]
    destroyed_sat_rs = three_rs[dist_mask]

    #grab the top ten (will there even be ten after all the cuts?)
    top_ten_mask = np.argsort(destroyed_sat_three_mass)[::-1][:10]
    top_ten_mass = destroyed_sat_three_mass[top_ten_mask]
    top_ten_x, top_ten_y, top_ten_z = destroyed_sat_three_x[top_ten_mask]-mmp_x, destroyed_sat_three_y[top_ten_mask]-mmp_y, destroyed_sat_three_z[top_ten_mask]-mmp_z
    top_ten_vx, top_ten_vy, top_ten_vz = destroyed_sat_three_vx[top_ten_mask]-mmp_vx, destroyed_sat_three_vy[top_ten_mask]-mmp_vy, destroyed_sat_three_vz[top_ten_mask]-mmp_vz
    top_ten_dist = destroyed_sat_dist[top_ten_mask]
    top_ten_rvir = destroyed_sat_rvir[top_ten_mask]
    top_ten_rs = destroyed_sat_rs[top_ten_mask]

    print(top_ten_dist)
    #Okay this works, now I just need to do a few tests and package it up

    f_out = np.zeros((10,9))
    f_out[:,0] = top_ten_mass
    f_out[:,1] = top_ten_x
    f_out[:,2] = top_ten_y
    f_out[:,3] = top_ten_z
    f_out[:,4] = top_ten_vx
    f_out[:,5] = top_ten_vy
    f_out[:,6] = top_ten_vz
    f_out[:,7] = top_ten_rvir
    f_out[:,8] = top_ten_rs
    
    np.savetxt('./'+output_dir+'/top_ten_'+str(halo_number)+'_'+sim_type+'.txt',f_out,header='# halo: {}, host mass at z = 3: {}*1.0e10 Msun, host Rvir at z = 3: {} kpc, host Rs at z = 3 {} \n#Mvir (Msun) x (kpc), y(kpc), z(kpc), vx (kms^-1), vy (kms^-1), vz (kms^-1), Rvir (kpc), Rs (kpc)'.format(halo_number,mmp_mass/1.0e10,mmp_rvir,mmp_rs))
