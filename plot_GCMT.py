#!/usr/bin/env python3

# This script loads all focal mechanisms from the GCMT catalogue (1976-2017) and 
# plots them

# %% codecell
from obspy.imaging.beachball import beach
from obspy import read_events, UTCDateTime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
#%matplotlib inline

# %% codecell
mindepth = 0; # [km] minimum depth
maxdepth = 900; # [km] maximum depth
minmag = 0; # minimum magnitude
maxmag = 10; # maximum magnitude
tstart = "2020-09-18T00:00:00" # start time 
tend =   "2020-09-25T00:00:00" # end time
latmin = -90 # minimum latitude of search region
latmax = 90 # maximum longitude of search region
lonmin = -180 # minimum longitude [-180, 180]
lonmax = 180 # maximum longitude [-180, 180]

figdir = './figs/' # directory for saving figures

# %% codecell

# Load quick CMT events from .ndk (https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk)
catCMT = read_events("https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk")

# # Load all earthquakes from .ndk (https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/jan76_dec17.ndk)
# catCMT = read_events('https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/jan76_dec17.ndk')

# Alternatively, can save the ndk file locally for quicker/offline access. In that case, uncomment following line
# catCMT = read_events('jan76_dec17.ndk')

# %% codecell
# Index earthquakes of interest
cat_filt = catCMT.filter("depth >= "+str(mindepth*1000),
                         "depth <= "+str(maxdepth*1000),
                         "magnitude >= "+str(minmag),
                         "magnitude <= "+str(maxmag),
                         "time >= "+str(UTCDateTime(tstart)),
                         "time <= "+str(UTCDateTime(tend)),
                         "latitude >= "+str(latmin),
                         "latitude <= "+str(latmax),
                         "longitude >= "+str(lonmin),
                         "longitude <= "+str(lonmax))
                             
# %% codecell
cat_filt
# print(catCMT_quick.__str__(print_all=True))
# cat_filt[3].origins[1].depth/1000

# %% codecell
# Plot focal mechanisms
    
llcrnrlat=latmin
urcrnrlat=latmax
llcrnrlon=lonmin
urcrnrlon=lonmax
parspace = 30 # spacing between lat lon labels

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
m = Basemap(projection='cyl', llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,
            llcrnrlon=0,urcrnrlon=360,
            resolution='c')

m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-90., 120., parspace),labels=[1,0,0,0],fontsize=12,linewidth=1)
m.drawmeridians(np.arange(0., 420., parspace),labels=[0,0,0,1],fontsize=12, linewidth=1)
m.drawmapboundary()
plt.title(tstart+' to '+tend,fontsize=18)

obj = {}
for iev,ev in enumerate(cat_filt):
    mt = ev.focal_mechanisms[0].moment_tensor.tensor
    mt_vec = [mt.m_rr, mt.m_tt, mt.m_pp, mt.m_rt, mt.m_rp, mt.m_tp ]
    lon = ev.origins[1].longitude
    lat = ev.origins[1].latitude
    if lon < 0:
        lon = lon + 360
    x, y = m(lon,lat)
    width = ((llcrnrlat-urcrnrlat)**2 +(llcrnrlon-urcrnrlon)**2)**0.5 / 15 * np.sqrt(ev.magnitudes[0].mag)/5 # width of beachballs in degrees
    obj[iev] = beach(mt_vec, xy=(x, y), width=width, linewidth=1, alpha=1, facecolor='r')
    obj[iev].set_zorder(10)
    ax.add_collection(obj[iev])
#     m.scatter(x,y,s=180,marker='o',color=(204/255, 0/255, 0/255),linewidths=1,edgecolors='k',zorder=10)
if not os.path.exists(figdir):
    os.makedirs(figdir)
fig.savefig(figdir+"GCMTevents.pdf", bbox_inches="tight")

