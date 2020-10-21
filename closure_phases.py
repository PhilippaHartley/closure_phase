from scipy.fftpack import ifftn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from math import sin, cos, sqrt, atan2, radians
from PIL import Image
from astropy.io import fits


mpl.rcParams['image.cmap'] = 'jet'



def mkgauss (naxes,pos,flux,fwhm,axrat=1.0,angle=0.0,ignore=4.0,dodist=False):

# note that total flux = peak flux in a pixel * 1.1331*FWHM**2
# angle is major axis East of North
    a = np.zeros (naxes[0]*naxes[1]).reshape(naxes[1],naxes[0])
    fwhm /= 1.66667
    if axrat==1.0 and angle==0.0:
        for i in range (naxes[1]):
            ydist=float(i)-pos[1]
            for j in range (naxes[0]):
                xdist=float(j)-pos[0]
                if xdist*xdist+ydist*ydist>ignore*ignore*fwhm*fwhm:
                    continue
                if not dodist:
                    a[i,j] = flux*np.exp(-(xdist*xdist+ydist*ydist)/ \
                                    (fwhm*fwhm))/(fwhm*fwhm*np.pi)
                else:
                    a[i,j] = np.hypot(xdist,ydist)
        return a
    sinth = np.sin(angle*np.pi/180.0)
    costh = np.cos(angle*np.pi/180.0)
    r = np.array([-sinth,costh,-costh,-sinth])
    rt = np.array([-sinth,-costh,costh,-sinth])
    sig = np.array([fwhm,0.0,0.0,fwhm*axrat])
    scr1 = mxmul (sig,r)
    scr2 = mxmul (rt, scr1)
    scr1 = mxinv (scr2)
    for i in range(naxes[1]):
        ydist=float(i)-pos[1]
        if abs(ydist)>ignore*fwhm:
            continue
        for j in range (naxes[0]):
            xdist = float(j) - pos[0]
            if abs(xdist)>ignore*fwhm:
                continue
            ex = scr1[0]*xdist+scr1[1]*ydist
            ey = scr1[2]*xdist+scr1[3]*ydist
            if not dodist:
                a[i,j] = (flux/axrat)*np.exp(-(ex*ex+ey*ey))/(fwhm*fwhm*np.pi)
            else:
                a[i,j] = np.hypot(ex,ey)/1.6666667

    return a 

def make_dynspec(lat1,lon1,lat2,lon2,Z):
    uv_centre = int(len(Z)/2)
    print (uv_centre)
    R = 6371
 
    raddec = dec*rad

    x1 = R*np.cos(lon1)*np.cos(lat1) 
    y1 = R*np.sin(lon1)*np.cos(lat1)
    z1 = R*np.sin(lat1)
    x2 = R*np.cos(lon2)*np.cos(lat2)
    y2 = R*np.sin(lon2)*np.cos(lat2)
    z2 = R*np.sin(lat2)
    dx = x2-x1
    dy = y2-y1
    dz = z2-z1
    us = np.array([])
    vs = np.array([])

    lamrang = np.linspace(c/sfreq,c/(sfreq+((channels+1)*ch_width)),channels)  

    print (lamrang)
    hrang =  np.linspace(0,(n_intervals+1)*time_int,n_intervals)
    print (hrang)
    dynspec = np.zeros((len(lamrang), len(hrang)))

    phase_plot = np.angle(Z)
    for lamcount, lam in enumerate(lamrang):
        
        # hour angle begins from same value each time, as if the sources all have same RA
        for hcount,h in enumerate(hrang):
            h = h*rad
            u = np.sin(h)*dx+np.cos(h)*dy
            v = -np.sin(raddec)*np.cos(h)*dx+np.sin(raddec)*np.sin(h)*dy+np.cos(raddec)*dz
            blength1 = (np.hypot(u,v))# in km!

       
            theta1rad = lam/(blength1*1000)
            theta1arcsec = theta1rad*rad2arcsec
           

            u*=1000 # in m
            v*=1000

            u/=lam # in lam
            v/=lam #      

          

            # each pixel of the fft2 is a frequemcy mode, with mode 0 at centre and nyquist mode at edge
            # maxfreq = 1/pixsize = eg 1/(0.01/rad2arcsec) # pix in radians

            umax = 1/(pix_size/rad2arcsec)
            vmax = 1/(pix_size/rad2arcsec)
            umax/=1000 
            vmax/=1000
           
            u = u/vmax
            v = v/vmax

            u = (u).astype(np.int)+uv_centre
            v = (v).astype(np.int)+uv_centre
          
            us = np.append(us,u)
            vs = np.append(vs,v)          
            phas = phase_plot[u,v]
       
            dynspec[lamcount,hcount]+=phas

    us = us.astype(np.int)
    vs = vs.astype(np.int)

    return dynspec, us, vs
    


def make_sources(x1,y1,x2,y2,dec,num):

    rad = (2.*np.pi)/360.
    xf = np.zeros((N, N))
    source_box_size = [101,101]
  #  g = mkgauss([source_box_size[0],source_box_size[1]], [(source_box_size[0]-1)/2-1,(source_box_size[1]-1)/2], 1, 3)
    g = mkgauss([101,101], [50,50], 1, 1)
    print (x1, x2)
    blc1 = [int(x1-source_box_size[0]/2),int(y1-source_box_size[1]/2)]
    trc1 = [int(x1+source_box_size[0]/2),int(y1+source_box_size[1]/2)]
    blc2 = [int(x2-source_box_size[0]/2),int(y2-source_box_size[1]/2)]
    trc2 = [int(x2+source_box_size[0]/2),int(y2+source_box_size[1]/2)]


    # add gaussians
    xf[blc1[0]:trc1[0],blc1[1]:trc1[1]]+=g
    xf[blc2[0]:trc2[0],blc2[1]:trc2[1]]+=g


    # add single pixels
  #  xf[int(x1), int(y1)] = 1
  #  xf[int(x2), int(y2)] = 1

    
    # shift *before* doing the ifft: removes phase jumps
    Z = np.fft.ifftn(np.fft.ifftshift(xf))
  

    if doplot:
        plt.imshow(np.angle(Z))
        plt.show()
    
    fig = plt.figure(figsize=(5,5))
    plt.imshow(xf, cmap = 'gray_r',interpolation = 'none') 
    plt.axis("off")
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    fig = plt.gcf()
    fig.savefig('training_images/B/'+np.str(num)+'.png', box_inches='tight', dpi=64)
    plt.show()
    plt.close()

    

    # lat and longs of some radio telescopes

    # European VLBI network
    Effelsberg = [50.5248, 6.8836]
    Yebes = [40.5336, 3.1112]

    # eMERLIN stations nb may not be exact
    Goonhilly =[50.048056, -5.181944]
    Cambridge = [52.167, 0.037056]
    Jodrell1 = [53.2369, -2.3075]
    Jodrell2 =  [53.2340, -2.3039]
    Darnhall = [53.164, -2.547]
    Knockin = [52.794, -2.993]
    Pickmere = [53.288, -2.462]
    Defford = [52.087346, -2.119104 ]

    #a source at dec 90 will trace the most circular ellipse in uv
    
    lat1 = (Cambridge[0])*rad
    lon1 = (Cambridge[1])*rad

    lat2 = (Knockin[0])*rad
    lon2 = (Knockin[1])*rad

    lat3 = (Goonhilly[0])*rad
    lon3 = (Goonhilly[1])*rad

    p1,us1,vs1 = make_dynspec(lat1,lon1,lat2,lon2,Z)
    p2,us2,vs2 = make_dynspec(lat2,lon2,lat3,lon3,Z) 
    p3,us3,vs3 = make_dynspec(lat1,lon1,lat3,lon3,Z)
    closphasdynspec = p1+p2-p3    

    if doplot:
        plt.subplot(221)
        plt.imshow(p1)
        plt.colorbar()
        plt.subplot(222)
        plt.imshow(p2)
        plt.colorbar()
        plt.subplot(223)        
        plt.imshow(p3)
        plt.colorbar()
        plt.subplot(224)   
        plt.imshow(closphasdynspec)
        plt.colorbar()
        plt.show()             




    if doplot:
        # plot uv coverage of the three antennae (making three baselines)
        plt.scatter(us1,vs1)
        plt.scatter(us2,vs2)
        plt.scatter(us3,vs3)
        plt.show()
    
    fig = plt.figure(figsize=(5,5))
    plt.imshow(closphasdynspec.T, interpolation = 'none')
    plt.xlabel('time')
    plt.ylabel('frequency')    
    plt.axis("off")
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0) 
    fig = plt.gcf()
    fig.savefig('training_images/A/'+np.str(num)+'.png', box_inches='tight', dpi=64)
    plt.close()


doplot = 1

# manage unwrapping
discont_val = 0.1
wrap_factor = 10

# approximate radius of earth in km
R = 6373.0    
    
# set up frequency resolution and bandwidth
channels = 64
sfreq = 1.4e9 # GHz
ch_width = 0.008e9 # GHz

# set up time resolution and number of integration times
time_range = 6 # hours
hours2deg = 360/24
hrang_max = time_range*hours2deg
n_intervals = 64 # number of integrations
time_int = hrang_max/n_intervals

# some constants
rad2arcsec = 206265.
c = 3e8
rad = (2.*np.pi)/360.  


# image dimensions: also determines dimensions of uv plane (fft2 array)
pix_size = 0.01 # in arcscec; this will depend on array and should be < half of max resolution of array (0.2arcsec for 1.4 GHz eMERLIN)
N = 1025

# set up source info
max_source_gap = 1 # arcsec
n_sources = 2


# select required number of image pairs
N_images = 1

# seed the random numbers to get reproducible results
np.random.seed(1)
for i in range(N_images):
    randoms = np.random.rand(4,N_images)
    # put sources in central  (2 arcsec)
    x1 = (N/2-(max_source_gap/2)/pix_size)+np.int(randoms[0,i]*(max_source_gap/pix_size))
    y1 = (N/2-(max_source_gap/2)/pix_size)+np.int(randoms[1,i]*(max_source_gap/pix_size))
    x2 = (N/2-(max_source_gap/2)/pix_size)+np.int(randoms[2,i]*(max_source_gap/pix_size))
    y2 = (N/2-(max_source_gap/2)/pix_size)+np.int(randoms[3,i]*(max_source_gap/pix_size))
    dec = 60

    # test with one source in centre and one located max_source_gap away
    x1 = (N/2)
    y1 = (N/2)
    x2 = (N/2)
    y2 = (N/2+(max_source_gap)/pix_size)

  
    make_sources(x1,y1,x2,y2,dec,i)
    
    


