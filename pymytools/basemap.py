from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from matplotlib.colors import LightSource, Normalize
from matplotlib.cm import ScalarMappable
from netCDF4 import Dataset
import pycpt
import obspy
import numpy as np
import configparser
import os

class myBaseMap(object):
    readConfig = False
    
    def _get_base(self, minlon, maxlon, minlat, maxlat, projection='lambert', resolution='i',**kwargs):
        """ Get basemap from Basemap of mpl_toolkits.basemap
        """
        lat_c = (maxlat+minlat)/2.0
        lon_c = (maxlon+minlon)/2.0
        if projection=='merc':
            m=Basemap(projection='merc', llcrnrlat=minlat-5., urcrnrlat=maxlat+5., llcrnrlon=minlon-5., \
                    urcrnrlon=maxlon+5., lat_ts=20, resolution=resolution, **kwargs)
        elif projection=='global':
            m=Basemap(projection='ortho',lon_0=lon_c, lat_0=lat_c, resolution=resolution, **kwargs)
        elif projection=='regional_ortho':
            m1 = Basemap(projection='ortho', lon_0=minlon, lat_0=minlat, resolution='l')
            m = Basemap(projection='ortho', lon_0=minlon, lat_0=minlat, resolution=resolution,\
                    llcrnrx=0., llcrnry=0., urcrnrx=m1.urcrnrx/mapfactor, urcrnry=m1.urcrnry/3.5,**kwargs)
        elif projection=='lambert':
            distEW, az, baz=obspy.geodetics.gps2dist_azimuth(minlat, minlon, minlat, maxlon+1) # distance is in m
            distNS, az, baz=obspy.geodetics.gps2dist_azimuth(minlat, minlon, maxlat+1., minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS, rsphere=(6378137.00,6356752.3142), resolution=resolution, projection='lcc', \
                lat_1=minlat-1, lat_2=maxlat+1, lon_0=lon_c, lat_0=lat_c, **kwargs)
        return m
    
    def __init__(self, minlon, maxlon, minlat, maxlat, projection='lambert', resolution='i', **kwargs):
        self.minlon = minlon
        self.maxlon = maxlon
        self.minlat = minlat
        self.maxlat = maxlat
        self.m      = self._get_base(minlon, maxlon, minlat, maxlat, projection, resolution, **kwargs)
        return

    def get_basemap(self, plate_kwargs={'draw': True}, etopo_kwargs={'draw': False},
                    lines_kwargs={'draw':True}, **kwargs):
        """ Get basemap for future plotting
        Parameters:
        ----------------------------------------------------------------------------------
        plate_kwargs -- parameters for self.draw_shp() (dictionary)
        etopo_kwargs -- parameters for self.draw_etopo() (dictionary)
        lines_kwargs -- parameters for self.draw_lines() (dictionary)
        dlat         -- latitude spacing in grid lines, default=2.
        dlon         -- longitude spacing in grid lines, default=2.
        latlabels    -- location for latitude ticks, default=[1,0,0,0] (W)
        lonlabels    -- location for longitude ticks, default=[0,0,1,0] (N)
        dashes       -- dash pattern for grids (default [1,1], i.e. 1 pixel on, 1 pixel off).
        cordft       -- coordinates label fontsize
        cordlw       -- grid line width
        coastlw      -- coast line width
        countrylw    -- country border line width
        statelw      -- state border line width
        ----------------------------------------------------------------------------------
        """
        dlat      = kwargs.get('dlat', 2.)
        dlon      = kwargs.get('dlon', 2.) 
        latlabels = kwargs.get('latlabels', [1,0,0,0])
        lonlabels = kwargs.get('lonlabels', [0,0,1,0])
        dashes    = kwargs.get('dashes', [2,2]) # dash pattern for meridians (default [1,1], i.e. 1 pixel on, 1 pixel off).
        cordft    = kwargs.get('cordft', 10) # coordinates label fontsize
        cordlw    = kwargs.get('cordlw', 1.0)
        coastlw   = kwargs.get('coastlw', 1.0)
        countrylw = kwargs.get('countrylw', 1.0)
        statelw   = kwargs.get('statelw', 1.0)
        zorder    = kwargs.get('zorder', 99)
        self.m.drawparallels(np.arange(-80.0, 80.0, dlat), linewidth=cordlw, dashes=dashes, labels=latlabels, fontsize=cordft, zorder=zorder)
        self.m.drawmeridians(np.arange(-180.0, 180.0, dlon), linewidth=cordlw, dashes=dashes, labels=lonlabels, fontsize=cordft, zorder=zorder)
        self.m.drawcoastlines(linewidth=coastlw, zorder=zorder)
        self.m.drawcountries(linewidth=countrylw, zorder=zorder)
        self.m.drawstates(linewidth=statelw)
        if plate_kwargs.get('draw'): # draw plate boundary
            self.draw_shp(**plate_kwargs)
        if lines_kwargs.get('draw'):
            self.draw_lines(**lines_kwargs)
        if etopo_kwargs.get('draw'):
            self.draw_etopo(**etopo_kwargs)
        pass

    def draw_shp(self, **kwargs):
        """ Draw shape file
        Prameters:
        ----------------------------------------------------------------------------------
        shpfname -- shape file name
        lw       -- line width, default orange
        color    -- line color, default 1.
        ----------------------------------------------------------------------------------
        """
        kwargs.pop('draw', None)
        shpfname = kwargs.pop("shpfname", '/home/wang/Code/ToolKit/Models/Plates/PB2002_boundaries')
        lw       = kwargs.pop("lw", 1.)
        color    = kwargs.pop("color", 'orange')
        try:
            self.m.readshapefile(shpfname, name=os.path.basename(shpfname), drawbounds=True, linewidth=lw, color=color, **kwargs) # draw plate boundary on basemap
        except IOError:
            print(f"Shape file {shpfname} doesn't exist! Continue without drawing plateboundaries.")
        return

    def draw_etopo(self, cbar_kwargs={'draw': False}, **kwargs):
        """ Draw etopo dataset on Basemap
        Parameter:
        ----------------------------------------------------------------------------------
        cbar_kwargs -- parameters for self.addColorBar() (dictionary)
        etopof      -- file name for the etopo dataset
        etopoCPT    -- cpt file for colormap
        azdeg       -- azimuth in degree for LightSource
        altdeg      -- altitude in degree for LightSource
        ----------------------------------------------------------------------------------
        """
        kwargs.pop('draw', None)
        if cbar_kwargs.get('draw'):
            addCbar = True
        else:
            addCbar = False
        etopof  = kwargs.get('etopof', '/home/wang/Code/ToolKit/Models/ETOPO1/ETOPO1_Ice_g_gmt4.grd')
        etopoCPT = kwargs.get('cpt', '/home/wang/Code/ToolKit/Models/ETOPO1/ETOPO1.cpt')
        azdeg   = kwargs.get('azdeg', 315.)
        altdeg  = kwargs.get('altdeg', 45.)
        try:
            mycm=pycpt.load.gmtColormap(etopoCPT)
            etopo1 = Dataset(etopof, 'r') # read in the etopo1 file which was used as the basemap
        except IOError:
            print("Couldn't read etopo data or color map file! Check file directory!")
            return
        lons = etopo1.variables["x"][:]
        west = lons<0 # mask array with negetive longitudes
        west = 360.*west*np.ones(len(lons))
        lons = lons+west
        lats = etopo1.variables["y"][:]
        z = etopo1.variables["z"][:]
        etopoz=z[(lats>(self.minlat-2))*(lats<(self.maxlat+2)), :]
        etopoz=etopoz[:, (lons>(self.minlon-2))*(lons<(self.maxlon+2))]
        lats=lats[(lats>(self.minlat-2))*(lats<(self.maxlat+2))]
        lons=lons[(lons>(self.minlon-2))*(lons<(self.maxlon+2))]
        lons = lons-360*(lons>180)*np.ones(len(lons))
        ind_lats = np.argsort(lats)
        ind_lons = np.argsort(lons)
        etopo_sort = etopoz[:,ind_lons]
        etopo_sort = etopoz[ind_lats,:]
        etopoZ = self.m.transform_scalar(etopo_sort, lons[ind_lons], lats[ind_lats], etopoz.shape[0], etopoz.shape[1]) # tranform the altitude grid into the projected coordinate
        ls = LightSource(azdeg=azdeg, altdeg=altdeg)
        vmin = etopoZ.min()
        vmax = vmin*85./(-110.)
        rgb = ls.shade(etopoZ, cmap=mycm, vmin=vmin,vmax=vmax, vert_exag=0.05, blend_mode='soft')
        img = self.m.imshow(rgb)
        if addCbar:
            cbar_kwargs.update({'label': "Topography (km)"})
            cbar_kwargs.update({'scale': 1000.})
            cbar_kwargs.update({'dv': 1000.})
            cbar_kwargs.update({'vmin': vmin})
            cbar_kwargs.update({'vmax': vmax})
            self.addColorBar(mycm, norm=True, **cbar_kwargs)
        return
    
    def addColorBar(self, cmap, img=None, norm=True, **kwargs):
        """ Add colorbar to basemap
        Parameter:
        ----------------------------------------------------------------------------------
        cmap   -- color map
        img    -- image the color map based on, only useful when norm is False
        norm   -- normalized the colorbar or not
        loc    -- location for the colorbar, default='bottom'
        size   -- size of the colorbar, default='3%'
        pad    -- colorbar pad, default='2%'
        vmin   -- minimum value on colorbar
        vmax   -- minimum value on colorbar
        dv     -- distance of the ticks
        scale  -- tick value scaling
        label  -- colorbar label
        ----------------------------------------------------------------------------------
        """
        loc   = kwargs.get('loc', 'bottom')
        size  = kwargs.get('size', "3%")
        pad   = kwargs.get('pad', "2%")
        vmin  = kwargs.get('vmin', None)
        vmax  = kwargs.get('vmax', None)
        dv    = kwargs.get('dv', None)
        scale = kwargs.get('scale', 1.)
        label = kwargs.get('label', '')
        if norm:
            sm = ScalarMappable(norm=Normalize(vmin=vmin,vmax=vmax),cmap=cmap)
            sm.set_array([])
            cbar = self.m.colorbar(sm, loc, size=size, pad=pad, extend='both', extendfrac='auto', extendrect='True')
        else:
            cbar = self.m.colorbar(img, loc, size=size, pad=pad, extend='both', extendfrac='auto', extendrect='True')
        if vmin is not None and vmax is not None:
            ticks = np.arange(np.ceil(vmin/scale), np.floor(vmax/scale)+0.2*dv/scale, dv/scale)
            tlabs = ['{:g}'.format(x) for x in ticks]
            cbar.set_ticks(ticks*scale)
            cbar.set_ticklabels(tlabs)
        cbar.set_label(label, rotation=0)
        self.m.drawmapboundary(fill_color="white")
        return
        

    def draw_lines(self, linef='/home/wang/Code/ToolKit/Models/wus_province_II.dat2', **kwargs):
        """ Draw line segments on Basemap, Segments of lines are represented by longitude latitude pairs
        (when lonlat True else the other way around)
        Paramenters:
        ----------------------------------------------------------------------------------
        linef  -- file name contaning the lines information
        deli   -- delimeter for different segments of lines
        color  -- line color on map
        lw     -- line with on map
        lonlat -- if True lon lat pair, else lat lon pair in input file
        ----------------------------------------------------------------------------------
        """
        kwargs.pop('draw', None)
        deli   = kwargs.get('deli', '>')
        color  = kwargs.get('color', 'r')
        lw     = kwargs.get('lw', 1.)
        lonlat = kwargs.get('lonlat', True)
        with open(linef) as fin:
            lines = fin.readlines()
        lons = np.array([], dtype=np.float64)
        lats = np.array([], dtype=np.float64)
        for iline in lines:
            if deli in iline:
                self.m.plot(lons, lats, c=color, latlon=True, lw=lw)
                lons = np.array([], dtype=np.float64)
                lats = np.array([], dtype=np.float64)
            else:
                if lonlat:
                    ilon, ilat = np.array(iline.split()).astype(np.float)
                else:
                    ilat, ilon = np.array(iline.split()).astype(np.float)
                lons = np.append(lons, ilon)
                lats = np.append(lats, ilat)
        self.m.plot(lons, lats, c=color, latlon=True, lw=lw)
        return
    
    def plot_markers(self, inarr, lonlat=True, marker='^', color='g', ms=9):
        """ Plot markers at longitude latitude (when lonlat True else the other way around) pair locations
        Paramenters:
        ----------------------------------------------------------------------------------
        inarr  -- input longitude latitude array (or the other way around)
        lonlat -- if True lon lat pair, else lat lon pair in input file
        marker -- shape of the marker, default='^'
        color  -- 
        ms     -- size of the marker, default=9
        ----------------------------------------------------------------------------------
        """
        lons = inarr[:,0]
        lats = inarr[:,1]
        if not lonlat:
            lons, lats = lats, lons
        lons[lons<0.] += 360.
        xx, yy = self.m(lons, lats)
        self.scatter(xx, xy, marker=marker, color='g', s=ms)
