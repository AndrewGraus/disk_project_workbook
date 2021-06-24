import os
import numpy as np

#I'm doing something else here and grabbing the 
#from the merger tree files, there are less halos
#in the merger tree for a reason that I don't understand
#I'm hoping its because there very small halos are cut out

#I think I need to revert to using the tree.dat file 

a_target = 1.0/(1.0+3.0) #want to find halos near z=3

for data_file in os.listdir('./disk_file_full/'):
    halo_number = data_file.split('_')[1]

    print(data_file)
    
    #loading data
    f = np.loadtxt('./full_mtrees/'+data_file+'/tree_0_0_0.dat',skiprows=46)

    Mvir = f[:,10]
    Rvir = f[:,11]
    cens = f[:,17:19]
    Vel = f[:,20:22]
    Vmax = f[:,16]
    ids_from_halo_cat = f[:,30]
    #Tree_id = snap['TreeNumber'][:] #not sure what equivalent is in tree file
    Tree_root_id = f[:,29]
    scale = f[:,0]

    print(np.unique(scale))

    scale_select = np.unique(scale)[np.argmin(np.abs(np.unique(scale)-a_target))]
    print('Looking for targets near {}'.format(scale_select))

    print(np.min(Tree_root_id))

    #Grab everything at that redshift that doesn't exist at z=0
    destroyed_cat_three = (Tree_root_id==0.0)&(scale==scale_select)
    everything_zero = (scale==1.0)
    print('Total halos: {}, destroyed halos at z=3: {}'.format(len(Mvir),np.sum(destroyed_cat_three)))

    cens_destroyed_three = cens[destroyed_cat_three]
    Mvir_destroyed_three = Mvir[destroyed_cat_three]
    vel_destroyed_three = Vel[destroyed_cat_three]
    
    Mvir_zero = Mvir[everything_zero]
    Tree_root_id_zero = Tree_root_id[everything_zero]
    

    #now I need the most massive halo at z = 0
    host_index = np.argmax(Mvir_zero)
    host_tree_id = Tree_root_id_zero[host_index] #This is the id of the host tree
    print('host_id: {}'.format(host_tree_id))
    
    Host_tree_mask = (Tree_root_id==host_tree_id)&(scale==scale_select) #find everything in the host tree near z=3
    #now grab positions, velocities and masses within that tree

    scale_host = scale[Host_tree_mask]
    cen_host = cens[Host_tree_mask]
    vel_host = Vel[Host_tree_mask]
    Mvir_host = Mvir[Host_tree_mask]
    Rvir_host = Rvir[Host_tree_mask]

    #now I have everything at z=3 and the host at z=3
    
    halo_cen_diff_three = cens_destroyed_three - cen_host
    halo_vel_diff_three = vel_destroyed_three - vel_host
    halo_dist_three = np.linalg.norm(sat_cen_diff_three)

    sat_three_mask = (halo_dist_three<Rvir_host)

    sat_mass_three = Mvir_destroyed_three[sat_three_mask]
    sat_cen_three = cen_destroyed_three[sat_three_mask]
    sat_vel_three = vel_destroyed_three[sat_three_mask]

    sat_ids_ten_largest = np.argsort(sat_mass_three)[::-1][:10]
    
    sat_mass_ten = sat_mass_three[sat_ids_ten_largest]
    sat_cen_ten = sat_cen_three[sat_ids_ten_largest]
    sat_vel_ten = sat_vel_three[sat_ids_ten_largest]
    
    top_ten_array = np.hstack((sat_mass_ten,sat_cen_ten,sat_vel_ten))
    np.savetxt('./top_ten_{}.txt'.format(data_file),top_ten_array)
