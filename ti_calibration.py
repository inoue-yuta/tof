# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:21:50 2018

@author: Inoue
"""

import numpy as np
import cv2
import os
import matplotlib.pyplot as plt

'''    path rule \\calibration dataset main\\each distance\\measured date and time\\*.dat
'''

calibration_dataset = 'cabe'
distances = ['def','defx2','long','longlong']

rows = 240
cols = 320
stat = 'median'
#stat = 'average'

def IQconvert(calibration_dataset = calibration_dataset, distances = distances, rows = rows, cols = cols, stat = stat):
    '''convert and save files
    '''
    frequencies = []
    freq_dirs = {}
    path = os.path.join(calibration_dataset,str(distances[0]))
    dirlist = os.listdir(path)
    path = os.path.join(path,dirlist[0])
    dirlist = os.listdir(path)
    for p in xrange(len(dirlist)):
        strp = dirlist[p]
        if 'MHz' in strp:
            if '.0' in strp:# or '30.0' in strp: # DEBUG
                f = strp[:strp.index('MHz')]
                frequencies.append(float(f))
                freq_dirs[float(f)] = strp
    print freq_dirs
    print 'Collect Frequencies'
    for i in xrange(len(distances)):
        calibration_data = []
        path = os.path.join(calibration_dataset,str(distances[i]))
        dirlist = os.listdir(path)
        print path
        for j in xrange(len(frequencies)):
            print(str(i)+':'+str(frequencies[j])+':'+str(distances[i]))
            calibration_data.append(os.path.join(path,dirlist[0],freq_dirs[frequencies[j]]))
            print 'Set path'
            fi = open(os.path.join(calibration_data[j], "ToFFrame_rawI_"+stat+".dat"), "rb")
            ibuf = fi.read(rows*cols*2)
            fq = open(os.path.join(calibration_data[j], "ToFFrame_rawQ_"+stat+".dat"), "rb")
            qbuf = fq.read(rows*cols*2)
            print 'Read file'
            rawIframe = np.frombuffer(bytearray(ibuf), np.int16, rows * cols)
            rawQframe = np.frombuffer(bytearray(qbuf), np.int16, rows * cols)        
            anglesRaw = np.arctan2(rawQframe.reshape(240, 320), rawIframe.reshape(240, 320))
            phase = (anglesRaw + (2 * np.pi)) % (2 * np.pi)
            amp = np.hypot(rawQframe.reshape(240, 320), rawIframe.reshape(240, 320))
            np.save(os.path.join(path, 'phase_'+stat), phase)
            np.save(os.path.join(path, 'amplitude_'+stat), amp)
            print 'Save .npy'
            cv2.imwrite(os.path.join(path, 'phase_'+stat+'_'+str(frequencies[j])+'.png'), np.uint8(phase * 255 / np.pi / 2.))
            cv2.imwrite(os.path.join(path, 'amplitude_'+stat+'_'+str(frequencies[j])+'.png'), np.uint8(amp * 255 / np.amax(amp)))
            np.savetxt(os.path.join(path, 'amplitude_'+stat+'_'+str(frequencies[j])+'.csv'),amp,delimiter=',')
            np.savetxt(os.path.join(path, 'phase_'+stat+'_'+str(frequencies[j])+'.csv'),phase,delimiter=',')
            depth = (300/frequencies[j])*(phase/(2.0*np.pi))/2.0
            np.savetxt(os.path.join(path, 'depth_'+stat+'_'+str(frequencies[j])+'.csv'),depth,delimiter=',')
            np.save(os.path.join(path, 'depth_'+stat+'_'+str(frequencies[j])+'.npy'),depth)
            print 'Convert depth'
    return

def Dcalibration(frequency = '80.0', stat = 'median'):
    depth_map = []
    for i in xrange(len(distances)):
        depth_map.append(np.load(os.path.join(calibration_dataset,str(distances[i]),'depth_'+stat+'_'+frequency+'.npy')).reshape(rows*cols))
    # graph image
    plt.figure(figsize=(8,6))
    plt.grid()
    xline = np.array(distances)
    print np.float64(xline)
    graph_map = np.array(depth_map).T
    print graph_map.shape
    param = np.zeros((2,rows*cols))
    for j in xrange(rows*cols):
        #plt.plot(np.float64(xline),graph_map[j])
        param[0,j], param[1,j] = np.polyfit(graph_map[j],np.float64(xline),1)
        #plt.plot(np.float64(xline),graph_map[j]*param[0,j]+param[1,j])
    #plt.legend(loc='best')
    #plt.title('sphere shape',fontsize = 24)    
    plt.xlabel("Frequency[MHz]",fontsize = 20)
    plt.ylabel("Depth[m]",fontsize = 20)
    plt.savefig(calibration_dataset+'\\depth_freq.png',format = 'png',dpi = 300)
    np.savetxt(calibration_dataset+'\\param_alpha_'+'.csv',np.float64(param[0].reshape(rows,cols)),fmt='%.5f',delimiter=',')
    np.save(calibration_dataset+'\\param_alpha_'+stat+'_'+'.npy',np.float64(param[0].reshape(rows,cols)))
    np.savetxt(calibration_dataset+'\\param_offset_'+'.csv',np.float64(param[1].reshape(rows,cols)),fmt='%.5f',delimiter=',')
    np.save(calibration_dataset+'\\param_offset_'+stat+'_'+'.npy',np.float64(param[1].reshape(rows,cols)))
    return

def Dadapting(frequency = '80.0', stat = 'median', mask = False):
    data = []
    mask_d = np.ones((rows,cols))*100
    if mask != False:
        mask_d = cv2.imread(calibration_dataset+'\\mask.png',0)
    alpha = np.load(os.path.join(calibration_dataset,'param_alpha_'+stat+'_.npy'))
    offset = np.load(os.path.join(calibration_dataset,'param_offset_'+stat+'_.npy'))
    for i in xrange(len(distances)):
        data.append(np.load(os.path.join(calibration_dataset,str(distances[i]),'depth_'+stat+'_'+frequency+'.npy')))
        data[i] = data[i]*alpha+offset
        data[i][mask_d<30] = 0
        np.savetxt(os.path.join(calibration_dataset,str(distances[i]),'depth_'+stat+'_'+frequency+'_calib.csv'),data[i],fmt='%.5f',delimiter=',')
        np.save(os.path.join(calibration_dataset,str(distances[i]),'depth_'+stat+'_'+frequency+'_calib.npy'),data[i])
        cv2.imwrite(os.path.join(calibration_dataset,str(distances[i]),'depth_'+stat+'_'+frequency+'_calib.png'),np.uint8(data[i] * 255 / np.max(data[i])))
    return data

def APread(num_distance = 8):
    '''read amp and depth files
    '''
    
    for i in xrange(len(distances)):
        print(str(i)+':'+str(distances[i]))
        path = os.path.join(calibration_dataset,str(distances[i]))
        dirlist = os.listdir(path)
        np.load()
    return amplitude_data, phase_data, depth_data

if __name__ == "__main__":
    IQconvert()
    #Dcalibration()
    Dadapting()
    #APread(8)
    