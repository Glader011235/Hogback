from __future__ import division
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:38:09 2017

@author: RCGlade
"""
"""
This code allows a hogback to evolve through time. A resistant layer of
rock, which weathers slowly, overlies a softer layer of rock that weathers
quickly. Resistant rock produces "blocks" which land on the adjoining
hillslope. Boundaries incise at a specified rate. User can set hogback
layer thickness, block size, and dip, as well as relative weathering and
incision rates. Trackable metrics included are time and space-averaged
slope, block height, weathering rate, and erosion rate. Parameter space
exploration allows use of different weathering rules and block movement
rules. Parameters that users need to specify are surrounded by many 
comment signs. 

R.C. Glade and R.S. Anderson
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

#################################################
#Initialize

#Hogback characteristics
thickness = 10
blockheight = 1
slopedegrees = 30
drop = 3
blockmove = 2

#weathering parameters
exponential = 1
wstar = 0.1
wdotnot = 1e-3
wdotnot_hard = wdotnot/100
H0 = 0

#for humped weathering curve
A = 0.1e-3
A_l1 = A/100
A_l2 = A/100

#transport parameters
k = 0.5
hstar = 0.2
edot_channel = 1e-4
Sc = 0.75
rockdensity = 2.3
soildensity = 2
stdev = 60

#time conditions
dt = 2
tmax = 1000000
t = np.arange(0,tmax+dt,dt)
imax = len(t)
nplots = 200
tplot = tmax/nplots
tchannel_stabilize = tmax
#################################################


#block release rules
hogbackslope = np.tan(np.deg2rad(slopedegrees))
dx = blockheight*np.cos(np.deg2rad(slopedegrees))
elevhogdrop = blockheight*np.sin(np.deg2rad(slopedegrees))
xchange = blockheight*np.cos(np.deg2rad(slopedegrees))

#domain setup
xmin = -250
xmax = 100
x = np.arange(xmin,xmax+dx,dx)
xedge = x[0:-1] + (dx/2)
zmax = 300
#xlimfixed = [min(x) max(x)]
#ylimfixed = [123 321]

#initial topography
fractureplanedegrees = 90 - slopedegrees
fractureplane = np.tan(np.deg2rad(fractureplanedegrees))
z1 = hogbackslope*x + zmax
z2 = -fractureplane*x + zmax
z2end = zmax-thickness*np.cos(np.deg2rad(slopedegrees))
z3 = np.ones(len(x))*z2end
zb = np.minimum(z1,z2)
zb = np.maximum(zb,z3)
ztopo0 = zb

#topography and lines for plotting
#intercepts
b1 = zmax
b2 = zmax - (thickness/(np.cos(np.deg2rad((slopedegrees)))))

#equations of lines bounding layers
z1l = (hogbackslope*x) + b1
z2l = z1 - (thickness/np.cos(np.deg2rad(slopedegrees)))
l2 = np.where((zb <= z1l) & (zb >= z2l))[0]

#block rules
blockend = len(x)-1
nrocks = thickness/blockheight
blockstart = l2[-1]+1
wdotdist = np.zeros(shape = 0,dtype=np.int64)

#initial soil conditions
H = np.zeros(len(x))
H[:] = H0

z = zb+H

#other initializations
#F(imax)=struct('cdata',[],'colormap',[])
wdottracking = np.empty(shape = 0)
peak = np.where(z == max(z[l2]))[0][0]
elevblock = np.zeros(len(wdotdist))
elevnext = np.zeros(len(wdotdist))
Hsave = np.empty(shape = 0)

#nframe = 0

#plot initial topography
fig1 = plt.figure(figsize = (6,4))
init_topo = plt.subplot()
plt.ion()
plt.show()
init_topo.plot(x,zb)
init_topo.axis('equal')

#initialize metric arrays
blockstart = 0
l = np.arange(blockstart+11,blockstart+81+5,5)
meanslope = np.zeros(len(l))
meanspacing = np.zeros(len(l))
meanH = np.zeros(len(l))
meanheight = np.zeros(len(l))
meanwdot = np.zeros(len(l))
meanflux = np.zeros(len(l))
erosion = np.zeros(len(l))
meanwdot2 = np.zeros(len(l))
analyticalslope = np.zeros(len(l))

meanslopesave = np.empty(shape = 0)
analyticalslopesave = np.empty(shape = 0)
meanfluxsave = np.empty(shape = 0)
meanHsave = np.empty(shape = 0)
meanheightsave = np.empty(shape = 0)
meanwdotsave = np.empty(shape = 0)
meanwdot2save = np.empty(shape = 0)
blocklocationsave = np.empty(shape = 0)
blockweathersave = np.empty(shape = 0)

#######################################

#run

for i in range(0,imax):
    
    #reevaluate location of hard layer
    l2 = np.where((zb <= z1l) & (zb >= z2l))[0]
    
    #set weathering rules
    #if using humped weathering rule
    if exponential == 0:
        
        wdot = wdotnot*np.exp(-H/wstar) + ((A*H/wstar))*np.exp(-H/wstar)
        wdot[l2] = 0
        
    
    
    #if using exponential weathering rule
    if exponential == 1:
        
        wdot = wdotnot * np.exp(-H/wstar)
        wdot[l2] = 0
        
  

    #moving blocks
    for j in range(0,len(wdotdist)):
        
        elevblock = z[wdotdist[j]]
        elevnext = z[wdotdist[j]+1]
        blockdrop = elevblock - elevnext
        
        if blockdrop >= 2*(blockheight-wdottracking[j]):
            
            wdotdist[j] = wdotdist[j] +1
            H[wdotdist[j]-1] = Hsave[j]
            zb[wdotdist[j]] = zb[wdotdist[j]]+(blockheight-(wdottracking[j]))+H[wdotdist[j]]
            zb[wdotdist[j]-1] = zb[wdotdist[j]-1] - (blockheight-(wdottracking[j])) - Hsave[j]
            Hsave[j] = H[wdotdist[j]]
            H[wdotdist[j]] = 0
            z = H+zb
                     
    #supply of blocks from hogback
    #keep track of elevation drop next to hogback
    elevdrop = z[l2[-1]] - z[l2[-1]+1]
    overshoot = elevdrop - drop

    #print '2 last zb = ', zb[-1]
    #if overshoot >= 2*drop:
        
        #disp('DANGER')
        
    peak = np.where(z == max(z[l2]))[0][0]

    #block release

    if ((elevdrop >= drop) & (t[i] > 0)):
        
        zmax = z[peak] - elevhogdrop
        zl2end = z[l2[-1]] - elevhogdrop
        zb[(peak-xchange/dx):(l2[-1]-xchange/dx)] = np.linspace(zmax,zl2end,((l2[-1]-xchange/dx)-(peak-xchange/dx)))
        zb[l2[-1]] = zb[l2[-1]] - drop
        l2 = np.where((zb >= z2l) & (zb <= z1l))[0]
        peak = np.where(z == max(z[l2]))[0][0]
        
        blockstart = l2[-1] +1
        
        #choose block distribution pattern
        #wdotdistnew=round(abs(stdev.*randn(1,nrocks))+blockstart); %use random distribution of block placement (right side of bell curve)
        #wdotdistnew=[blockstart:2:blockstart+((2*nrocks)-1)]; %place blocks immediately downslope, with one space in between each other
        wdotdistnew = np.arange(blockstart,blockstart+nrocks,1,dtype = np.int64)
        for_deletion = np.where(wdotdistnew>blockend)[0]
        wdotdistnew = np.delete(wdotdistnew, for_deletion, axis = 0)
        
        zb[wdotdistnew] = zb[wdotdistnew] + blockheight + H[wdotdistnew]
        Hsavenew = H[wdotdistnew]
        Hsave = np.concatenate((Hsave,Hsavenew))
        
        H[wdotdistnew] = 0
        wdotdist = np.concatenate((wdotdist,wdotdistnew),axis = 0)
        wdottracking = np.concatenate((wdottracking,np.zeros(len(wdotdistnew))),axis = 0)
     
    #print '3 last zb = ', zb[-1]    
    #update weathering array to keep track of block heights
    if exponential == 0:
        
        #humped weathering rule
        wdot_Hsave = wdotnot*np.exp((-Hsave+(blockheight-wdottracking))/wstar) + (((A_l1*Hsave)/wstar)*np.exp(-Hsave/wstar))
        Hsave = Hsave + (wdot_Hsave*dt)
        wdot[wdotdist] = wdotnot_hard*np.exp(-H[wdotdist]/wstar)+(((A_l2*H[wdotdist])/wstar)*np.exp(-H[wdotdist]/wstar))
        wdottracking = wdottracking + wdot[wdotdist]*dt
        
    if exponential == 1:
        
        #exponential weathering rule
        wdot_Hsave = wdotnot*np.exp(-(Hsave+(blockheight-wdottracking))/wstar)
        Hsave = Hsave+(wdot_Hsave*dt)
        wdot[wdotdist] = wdotnot_hard*np.exp(-H[wdotdist]/wstar)
        #wdot[wdotdist] = wdotnot_hard
        wdottracking = wdottracking + wdot[wdotdist]*dt
        
    for_soil_deletion = np.where((wdotdist == len(x)-2) | (wdotdist == len(x)-1))
    Hsave = np.delete(Hsave,for_soil_deletion,axis = 0)
    wdotdist = np.delete(wdotdist,for_soil_deletion,axis = 0)
    wdottracking = np.delete(wdottracking,for_soil_deletion,axis = 0)
    
    #soil transport
    slope = np.diff(z)/dx
    
    #find peaks and troughs to find H at usplope cell
    signdiff = np.sign(slope)
    xedge = x[0:-1]+(dx/2)
    lt = np.where(signdiff < 0)[0]
    rt = np.where(signdiff > 0)[0]
    
    #find soil thickness at upslope cell
    hedge = np.zeros(len(xedge))
    hedge[1] = H[2]
    hedge[-1] = H[-2]
    hedge[rt] = H[rt+1]
    hedge[lt] = H[lt]
    
    #flux from Johnstone and Hilley
    q = -k*slope*hstar*(1-(np.exp(-hedge/hstar)))
    
    #flux with nonlinear component
    #q = k*slope*hstar*(1-(np.exp(-hedge/hstar)))/(1-((np.abs(slope)/Sc)^2))
    
    q = np.insert(q,0,q[0])
    q = np.append(q,q[-1])
        
    #now conserve mass
    dh_dt = wdot - np.diff(q)/dx
    H = H+(dh_dt*dt)
    H = np.maximum(0,H)
    zb[1:-1] = zb[1:-1] - (wdot[1:-1]*dt)

    #print '4 last zb = ', zb[-1]
    #Get rid of blocks once they have completely weathered
    H[wdotdist[wdottracking >= blockheight]] = Hsave[wdottracking >= blockheight]
    zb[wdotdist[wdottracking >= blockheight]] = z[wdotdist[wdottracking >= \
        blockheight]] - H[wdotdist[wdottracking >= blockheight]]
    for_Hsave_deletion = np.where(wdottracking >= blockheight)[0]
    Hsave = np.delete(Hsave,for_Hsave_deletion,axis = 0)
    wdotdist = np.delete(wdotdist,for_Hsave_deletion,axis = 0)
    wdottracking = np.delete(wdottracking,for_Hsave_deletion, axis = 0)
    
    #now impose channel boundary condition
    if t[i] <= tchannel_stabilize:
        zb[0] = zb[0] - (edot_channel*dt)
        zb[-1] = zb[-1]-(edot_channel*dt)
        H[0] = 1
        H[-1] = 1
    else: #no further incision
        H[0] = 1
        H[-1] = 1
    zb_test = zb[-1]
  
    #print 'final last zb = ', zb[-1]    
    z = H+zb
    #plt.plot(x,z)
        
    #Metrics
    height = blockheight-wdottracking
    wdotdistsort = np.sort(wdotdist)
    wdotmean = wdot
    
	
    
    
    
    
#    if t[i] > 2000000:
#        
#        l = np.arange(blockstart+11,blockstart+81+5,5)
#        
#        for j in range(0,len(l)):
#            
#            meanslope[j] = ((z[l[j]] - ((np.sum(height[wdotdist==l[j]]))))-\
#                (z[l[j]+10]-(np.sum(height[wdotdist == l[j]+10]))))/\
#                (x[l[j]+10]-x[l[j]])
#             
#            meanspacing[j] = np.mean(np.diff((wdotdistsort[(wdotdistsort >= \
#                l[j]) & (wdotdistsort <= l[j]+10)])))            
#             
#            meanheight[j] = np.mean(height[(wdotdist >= l[j]) & (wdotdist <= \
#                l[j] +10)] )
#            
#            meanH[j] = np.mean(H[l[j]:l[j]+10])
#            
#            meanwdot[j] = np.mean(wdot[l[j]:l[j]+10])
#            
#            #meanwdot2(j) = nanmean(wdotmean(l(j):l(j)+10))
#            
#            meanflux[j] = np.mean(q[l[j]:l[j]+10])
#            
#            analyticalslope = meanheight/meanspacing + ((meanwdot*\
#                (np.abs((np.abs(np.arange(x[blockstart+16],x[blockstart+86+5],5))\
#                -np.abs(x[blockstart])))))/(k*hstar*(1-np.exp(-meanH/hstar))));
#                


