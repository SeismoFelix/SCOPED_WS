#Code for Running Specfem3D for the 6 elementary sources and then, processing the output files for being read by MTUQ. 
#Felix Rodriguez Cardozo, September 2022
#Please refer to the manual before using this algorithm. 

import argparse
from obspy.core import read
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
import glob
import os
import numpy as np

def check_solver(id,depth):
    
    if os.path.exists('SOLVER_REPO/{}/{}'.format(id,depth)):
        list_MT = glob.glob('SOLVER_REPO/{}/{}/M*'.format(id,depth))
        if len(list_MT) == 6:
            return(1)
        else:
            print('Elementary sources directories in SOLVER_REPO/{}/{}/ are not complete. You need 6: Mpp,Mrp,Mrr,Mrt,Mtp,Mtt'.format(id,depth))
            return(0)
    else:
        print('I cannot find the path SOLVER_REPO/{}/{}'.format(id,depth))
        return(0)

def check_synthetics(id,depth,MT):
    if os.path.exists('SOLVER_REPO/{}/eventinfo.txt'.format(id)):
        if os.path.exists('SOLVER_REPO/{}/STATIONS'.format(id)):
            cont = 0
            for m in MT:
                check_syn = glob.glob('SOLVER_REPO/{}/{}/{}/OUTPUT_FILES/*.sac'.format(id,depth,m))
                if len(check_syn) == 0:
                    print(' I cannot find synthetics in SOLVER_REPO/{}/{}/{}/OUTPUT_FILES/'.format(id,depth,m))
                    print('It seems you have not run the solver for {} elementary source'.format(m))
                else:
                    cont = cont+1
            if cont == 6:
                return(1)
            else:
                print('Your elementary source simulations are not complete')
                return(0)
        else:
            print('I cannot find the file SOLVER_REPO/{}/STATIONS'.format(id))
            print('Put there the STATIONS file used in Specmfem3D in Geographical Coordinates')
            return(0)
    else:
        print('I cannot find the  file SOLVER_REPO/{}/eventinfo.txt'.format(id))
        print('Make a file like this: YEAR-MONTH-DAYTHOUR:MINUTE:SECONDS LATITUDE LONGITUDE DEPTH MAGNITUDE.\ For example:')
        print('2017-12-01T02:32:44.0 30.7340 57.3900 15.0 5.9')
        return(0)

def process_elementary_sources(MT,id,depth):
    os.chdir('SOLVER_REPO/{}/{}'.format(id,depth))
    #print('###### Converting Specfem3D files into SAC and mini-seed files #####')
    mkdir_output = 'mkdir PROCESSED'
    os.system(mkdir_output)

    for m in MT:
        print('**********Renaming SAC files and getting the scale factor from the CMTSOLUTION file for {}********'.format(m))
        sc = rename_get_sc(m)
        print('********Rotating to R and T {} SAC files ****'.format(m))    
        rotate(m)

    print('*******Scaling Green Functions*****') 
    scale_gf(sc)
    print('*******Resampling Green Functions*****')
    new_delta = 0.2
    change_sampling_rate(new_delta)

def rename_get_sc(m):

    DIR = '{}/OUTPUT_FILES/*.sac'.format(m)

    list_sac=glob.glob(DIR)

    #Read CMTSOLUTION file for guessing the scaling factor
    open_cmtsolution= open('{}/DATA/CMTSOLUTION'.format(m),'r')
    read_cmtsolution = open_cmtsolution.readlines()
    
    for line in read_cmtsolution:
        aux = (line.split()[0].split(':')[0]).upper()
        #print(aux,m)
        if aux == m:
            sa = float(line.split()[-1])
            sc = (sa/1e+7)

    for i in range(len(list_sac)):
        comp_name = list_sac[i].split('/')[-1].split('.')[2][-1]
        #print(comp_name)
        st = read(list_sac[i])
        #The network correction is just because I made a mistake in the STATION file in Specfem3D 
        #st[0].stats.sac.knetwk = 'IR'
        #st[0].stats.network = 'IR'
        st[0].stats.sac.khole = ''
        st[0].stats.location = ''
        #st[0].stats.sac.kcmpnm = 'BH{}'.format(comp_name)
        #st[0].stats.channel = 'BH{}'.format(comp_name)

        new_name = rename(list_sac[i],m)
        st.write(new_name, format='SAC')
        #print (st)
        mv_sac = 'mv {} PROCESSED/'.format(new_name)
        #print(mv_sac)
        os.system(mv_sac)
        #print('\n')

    return(sc)

def rename(old_name,m):
    #Mrp/OUTPUT_FILES/IC.VRN.BXZ.sem.sac
    #Mrp/OUTPUT_FILES/IR.VRN..Z.Mrp.sac
    
    aux=old_name.split('/')[-1].split('.')
    #network = aux[0]
    network = aux[0]
    st = aux[1]
    comp = aux[2][-1]
    if comp == 'X':
        comp = 'E'
    elif comp == 'Y':
        comp = 'N'
    new_name = '{}/OUTPUT_FILES/{}.{}..{}.M{}.sac'.format(m,network,st,comp,m[1:3].lower())
    return(new_name)

def rotate(m):
    #PROCESSED/IR.CHBR..Z.Mrp.sac
    list_st = glob.glob('PROCESSED/*..Z.*.sac'.format(m))
    for t in list_st:
        stZ = read(t)
        net = stZ[0].stats.network
        source_lat = stZ[0].stats.sac.stla
        source_lon = stZ[0].stats.sac.stlo
        st_lat = stZ[0].stats.sac.evla
        st_lon = stZ[0].stats.sac.evlo
        baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)
        comp1 = stZ[0].stats.sac['kcmpnm'][0:-1]
        
        st_name = t.split('.')[1]
	    #print('Rotating for {},{}'.format(m,st_name))
        stN = read('PROCESSED/{}.{}..N.M{}.sac'.format(net,st_name,m[1:3].lower()))
        stE = read('PROCESSED/{}.{}..E.M{}.sac'.format(net,st_name,m[1:3].lower()))

        rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])
    
        stN[0].data = rotated[0]
        stN[0].stats.sac.kcmpnm = '{}R'.format(comp1)
        stN[0].stats.channel = '{}R'.format(comp1)
        stN[0].stats.sac.cmpaz = baz[2]

        stE[0].data = rotated[1]
        stE[0].stats.sac.kcmpnm = '{}T'.format(comp1)
        stE[0].stats.channel = '{}T'.format(comp1)

        if baz[2] + 90 >= 360:
            stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
        else :
            stE[0].stats.sac.cmpaz = baz[2] + 90

        r_name = 'PROCESSED/{}.{}..R.M{}.sac'.format(net,st_name,m[1:3].lower())
        t_name = 'PROCESSED/{}.{}..T.M{}.sac'.format(net,st_name,m[1:3].lower())
        
        stN.write(r_name,format='SAC')
        stE.write(t_name,format='SAC')

def scale_gf(sc):
    list_ev = glob.glob('PROCESSED/*.sac')

    for ev in list_ev:
        #print('Scaling {}'.format(ev))
        st = read(ev)
        st[0].data = st[0].data/sc
        st.write(ev,format='SAC')

def change_sampling_rate(new_delta):
    print('\tResampling the delta to be {}\n'.format(new_delta))
    data = glob.glob('PROCESSED/*.sac')
    print(len(data))

    for d in data:
        st = read(d)
        st.resample(sampling_rate=1/new_delta)
        st.write(d,format='SAC')

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-ev",type=str,help="Event ID. E.g. -ev  20171201023244")
    parser.add_argument("-ed",type=str,help="Earthquake depth in km E.g. -ed  6")

    args = parser.parse_args()
    id = str(args.ev)
    depth = str(args.ed)

    cwd = os.getcwd()
    MT = ['MPP','MRP','MRR','MRT','MTP','MTT']
    #MT = ['Mrp','Mrr','Mrt','Mtp','Mtt']
    go = check_solver(id,depth)

    if go == 0:
        print('End')
    else:

        print('Pre-processing Specfem3D output files')
        go2 = check_synthetics(id,depth,MT)
        if go2 == 0:
            print('End')
        else:
            process_elementary_sources(MT,id,depth)