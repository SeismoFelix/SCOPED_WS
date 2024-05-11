import argparse
from obspy.core import read
from obspy import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
import glob
import os
import numpy as np
from pandas import DataFrame

class Station:
    def __init__(self,name,lat,lon,network):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.network = network

def clean():
    print('rm -r PROCESSED')
    os.system('rm -r PROCESSED')
    
def get_event_info(path):
    print('\tGathering information about Specfem3D simulation\n')

    open_cmt = open('{}/DATA/CMTSOLUTION'.format(path),'r')
    header_cmt = open_cmt.readlines()[0].split()
    open_cmt.close()

    #print(header_cmt)
    year = int(header_cmt[1])

    month = int(header_cmt[2])
    if month < 10:
        month = '0{}'.format(month)

    day = int(header_cmt[3]) 
    if day < 10:
        day = '0{}'.format(day)

    hour = int(header_cmt[4])
    if hour < 10:
        hour= '0{}'.format(hour) 

    min = int(header_cmt[5])
    if min < 10:
        min= '0{}'.format(min) 

    sec = float(header_cmt[6])
    if sec < 10:
        sec= '0{}'.format(sec)

    ev_id = '{}{}{}{}{}{}'.format(year,month,day,hour,min,int(sec))

    #2017-12-01T02:32:44.0
    time = '{}-{}-{}T{}:{}:{}'.format(year,month,day,hour,min,sec)
    evla = float(header_cmt[7])
    evlo = float(header_cmt[8])
    evdp = float(header_cmt[9])

    return(evla,evlo,evdp,time,ev_id)

def convert_mtuq_format(path,ev_id,time):
    print('\tConverting plain text seismograms into SAC files\n')
    
    #evla,evlo,evdp,time,ev_id,utm_on = get_event_info(path)
    mkdir_output = 'mkdir PROCESSED'
    os.system(mkdir_output)

    DIR = '{}/OUTPUT_FILES/*.sac'.format(path)

    list_sac=glob.glob(DIR)


    #IR.GHIR.MXZ.sem.sac -> 20171201023244.IR.GHIR..MX.z

    for i in range(len(list_sac)):
        print('Renamig {}'.format(list_sac[i]))

        tr = read(list_sac[i])
        network = tr[0].stats.sac['knetwk']
        st = tr[0].stats.sac['kstnm']
        comp1 = tr[0].stats.sac['kcmpnm'][0:-1]
        comp2 = tr[0].stats.sac['kcmpnm'][-1].lower()
        new_name = 'PROCESSED/{}.{}.{}..{}.{}'.format(ev_id,network,st,comp1,comp2)
        tr[0].stats.sac.khole = ''
        tr[0].stats.location = ''

        tr.write(new_name, format='SAC')
        print(new_name)

def rotate(process_path,ev_id):
    print('\tRotating the radial and transverse\n')
    #20171201023244.IR.JHBN..HH.z
	
    list_st = glob.glob('{}/*.z'.format(process_path))
    for t in list_st:

        stZ = read(t)
        st_name = stZ[0].stats.station
        net = stZ[0].stats.network
        source_lat = stZ[0].stats.sac.stla
        source_lon = stZ[0].stats.sac.stlo
        st_lat = stZ[0].stats.sac.evla
        st_lon = stZ[0].stats.sac.evlo
        baz = gps2dist_azimuth(source_lat,source_lon,st_lat,st_lon)
        comp1 = stZ[0].stats.sac['kcmpnm'][0:-1]
        comp2 = stZ[0].stats.sac['kcmpnm'][-1].lower()
        
        
        print('Rotating for {},{}'.format(process_path,st_name))
        stN = read('{}/{}.{}.{}..{}.n'.format(process_path,ev_id,net,st_name,comp1))
        stE = read('{}/{}.{}.{}..{}.e'.format(process_path,ev_id,net,st_name,comp1))
        #print(stN)
        #print(stE)

        rotated = rotate_ne_rt(stN[0].data,stE[0].data,baz[1])
        
        stN[0].data = rotated[0]
        stN[0].stats.sac.kcmpnm = '{}R'.format(stN[0].stats.sac.kcmpnm[0:-1])
        stN[0].stats.channel = '{}R'.format(stN[0].stats.sac.kcmpnm[0:-1])
        stN[0].stats.sac.cmpaz = baz[2]

        stE[0].data = rotated[1]
        stE[0].stats.sac.kcmpnm = '{}T'.format(stE[0].stats.sac.kcmpnm[0:-1])
        stE[0].stats.channel = '{}T'.format(stE[0].stats.sac.kcmpnm[0:-1])

        if baz[2] + 90 >= 360:
            stE[0].stats.sac.cmpaz = baz[2] + 90 - 360
        else :
            stE[0].stats.sac.cmpaz = baz[2] + 90

        #20171201023244.IR.JHBN..HH.z
        r_name = '{}/{}.{}.{}..{}.{}'.format(process_path,ev_id,net,st_name,comp1,stN[0].stats.channel[-1].lower())
        t_name = '{}/{}.{}.{}..{}.{}'.format(process_path,ev_id,net,st_name,comp1,stE[0].stats.channel[-1].lower())
        print('\twriting {}'.format(r_name))
        print('\twriting {}'.format(t_name))
        #print('*****')
        #print(stN)
        #print(stE)
        stN[0].write(r_name,format='SAC')
        stE.write(t_name,format='SAC')

def padd_zeros(event,ev_id,time):
    print('\tAdding {}s of zeros to {} seismograms\n'.format(time,event))
    #id_event = event.split('/')[-1]
    data = glob.glob('{}/{}.*'.format(event,ev_id))
    extra_time = time

    for d in data:
        st = read(d)
        ns = round(extra_time/st[0].stats.delta)
        extra_samples = np.zeros(ns)
        st[0].data = np.concatenate((extra_samples,st[0].data,extra_samples))
        st[0].stats.starttime = st[0].stats.starttime - ns*st[0].stats.delta
        #st[0].data = 100*st[0].data 
        #st[0].stats.endtime = st[0].stats.endtime + extra_time
        st.write(d,format='SAC')

def scale_amplitude(event,scale,ev_id):
    print('\tMultiplying by {}, amplitude of {} seismograms\n'.format(scale,event))
    data = glob.glob('{}/{}*'.format(event,ev_id))
    print(len(data))

    for d in data:
        st = read(d)
        st[0].data = st[0].data*scale
        st.write(d,format='SAC')

def change_sampling_rate(event,ev_id,new_delta):
    print('\tResampling the delta to be {}\n'.format(new_delta))
    data = glob.glob('{}/{}*'.format(event,ev_id))
    print(len(data))

    for d in data:
        st = read(d)
        st.resample(sampling_rate=1/new_delta)
        st.write(d,format='SAC')


def write_weight_all(processed_dir,components):

    print('\tWriting weights.dat files for events in {} '.format(processed_dir))
    list_events = glob.glob('{}/*'.format(processed_dir))
     
    data = glob.glob('{}/*.r'.format(processed_dir))

    open_w_all = open('{}/weights.dat'.format(processed_dir),'w')

    for t in range(len(data)):
        name = data[t].split('/')[-1][0:-1]
        st = read(data[t])
        distance = gps2dist_azimuth(st[0].stats.sac['evla'],st[0].stats.sac['evlo'],st[0].stats.sac['stla'],st[0].stats.sac['stlo'])
        if t < len(data) - 1:
            line=' {} {} {}   0.0   0.0      0      0      0\n'.format(name,np.round(distance[0]/1000),components)
        else:
            line=' {} {} {}   0.0   0.0      0      0      0'.format(name,np.round(distance[0]/1000),components)
        open_w_all.write(line)
    open_w_all.close()

if __name__=='__main__':

    print('Read the notebook for more details')
   
    

    
