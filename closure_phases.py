from scipy.fftpack import ifftn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import sin, cos, sqrt, atan2, radians
from PIL import Image
from astropy.io import fits
def mkgauss2():

    x, y = np.meshgrid(np.linspace(-1,1,10), np.linspace(-1,1,10))
    print (x,y)
    d = np.sqrt(x*x+y*y)
    sigma, mu = 1.0, 0.5
    g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
    print("2D Gaussian-like array:")
    plt.imshow(g)
    plt.show()


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
    R = 6371
    rad = (2.*np.pi)/360.
    print (dec)
    raddec = dec*rad
    x1 = R*np.cos(lon1)*np.cos(lat1)
   # print R

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
    channels = 64
    sfreq = 1.4e9
    ch_width = 0.002e9
    lamrang = np.linspace(c/sfreq,c/(sfreq+((channels+1)*ch_width)),channels)  
    hrang =  np.arange(0,120,1)
    dynspec = np.zeros((len(lamrang), len(hrang)))
    print (dynspec.shape)
    phase_plot = np.angle(Z)
    for lamcount, lam in enumerate(lamrang):
        
        # hour angle begins from same value each time, as if the sources all have same RA
        for hcount,h in enumerate(hrang):
            h = h*rad
            u = np.sin(h)*dx+np.cos(h)*dy
            v = -np.sin(raddec)*np.cos(h)*dx+np.sin(raddec)*np.sin(h)*dy+np.cos(raddec)*dz
            blength1 = (np.hypot(u,v))# in km!!!
       
            theta1rad = lam/(blength1*1000)

            theta1arcsec = theta1rad*rad2arcsec
 

           
       

            u*=1000 # in m
            v*=1000

            u/=lam # in lam
            v/=lam #

            u/=1000# in klam
            v/=1000


            # each pixel of the fft2 is a frequemcy mode, with mode 0 at centre and nyquist mode at edge
            # maxfreq = 1/pixsize = 1/0.1 arcsec = 10?
            nyquist = 1/pix_size
          
            scale = nyquist/(len(Z)/2.)
            u*=scale
            v*=scale # deal with pix size?
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
    g = mkgauss([101,101], [50,50], 1, 5)


 

    blc1 = [int(x1-source_box_size[0]/2),int(y1-source_box_size[1]/2)]
    trc1 = [int(x1+source_box_size[0]/2),int(y1+source_box_size[1]/2)]
    blc2 = [int(x2-source_box_size[0]/2),int(y2-source_box_size[1]/2)]
    trc2 = [int(x2+source_box_size[0]/2),int(y2+source_box_size[1]/2)]



    xf[blc1[0]:trc1[0],blc1[1]:trc1[1]]+=g
    xf[blc2[0]:trc2[0],blc2[1]:trc2[1]]+=g
    #xf[x1, y1] = 10
    #xf[x1, y1] = 10



    # load real data
    #xf = fits.getdata('../closure_phases_casa/example.fits')



    if doplot:
        plt.imshow(xf)
        plt.show()
    Z = np.fft.fftshift(np.fft.fft2(xf))



    
    fig = plt.figure(figsize=(5,5))
    plt.imshow(xf, cmap = 'gray_r',interpolation = 'none') 
    plt.axis("off")
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.savefig('training_images/B/'+np.str(num)+'.png', box_inches='tight', dpi=64)
    plt.close()
    

    # lat and longs of some radio telescopes
    lat_longs = np.array([[53.244809, -2.309140],[52.920553, 6.603082],[50.529832, 6.882804],[49.145437, 12.878118],[53.090425, 18.558107],[57.396055, 11.926298],[31.104712, 121.198691],[40.5336, 3.1112],[43.820799, 87.618012],[60.218872, 24.393522],[60.535748, 29.781042],[43.831820, 41.588759],[51.772618, 102.238942],[57.549963, 22.130064],[39.503666, 9.247842],[36.883122, 14.986711],[-25.883272, 27.685024]])
    Effelsberg = [50.5248, 6.8836]
    Jodrell = [53.2369, -2.3075]
    Yebes = [40.5336, 3.1112]

    Goonhilly =[50.048056, -5.181944]
    Cambridge = [52.167, 0.037056]

    #a source at dec 90 will trace the most circular ellipse in uv



    # approximate radius of earth in km
    R = 6373.0
    
    lat1 = (Jodrell[0])*rad
    lon1 = (Jodrell[1])*rad
    lat3 = (Goonhilly[0])*rad
    lon3 = (Goonhilly[1])*rad
    lat2 = (Cambridge[0])*rad
    lon2 = (Cambridge[1])*rad
    
 

    p1,us1,vs1 = make_dynspec(lat1,lon1,lat2,lon2,Z)
    p2,us2,vs2 = make_dynspec(lat3,lon3,lat1,lon1,Z) 
    p3,us3,vs3 = make_dynspec(lat2,lon2,lat3,lon3,Z)


    plt.imshow(p1, cmap = 'jet')
    plt.show()
    plt.imshow(p2, cmap = 'jet')
    plt.show()
    plt.imshow(p3, cmap = 'jet')
    plt.show()


    wrap_factor = 10.
    print (p1.shape)
    for i in range(len(p1[:,0])):
        print (i)

        p1_i = np.unwrap(p1[i,:]*wrap_factor, discont = 0.1)/wrap_factor

        try:
            uw_p1_all =  np.vstack((uw_p1_all, p1_i))
        except:
            uw_p1_all =p1_i



    for i in range(len(uw_p1_all[0,:])):
        print (i)

        p1_i2 = np.unwrap(uw_p1_all[:,i]*wrap_factor, discont = 0.1)/wrap_factor

        try:
            uw_p1_all2 =  np.vstack((uw_p1_all2, p1_i2))
        except:
            uw_p1_all2 =p1_i2


    plt.clf()
    plt.imshow(uw_p1_all.T, origin = 'lower', cmap = 'jet')
    plt.colorbar()
 #   plt.show()
    #plt.savefig('unwrapped_dynspc_stacked.png')
    plt.clf() 

    plt.clf()
    plt.imshow(uw_p1_all2.T, origin = 'lower', cmap = 'jet')
    plt.colorbar()
    plt.show()
    #plt.savefig('unwrapped_dynspc_stacked.png')
    plt.clf() 


    '''
   
    closphasdynspec = p1+p2-p3 


    if doplot:
        # plot uv coverage of the three antennae (making three baselines)
        plt.scatter(us1,vs1)
        plt.scatter(us2,vs2)
        plt.scatter(us3,vs3)
        plt.show()
        plt.scatter(np.arange(len(closphasdynspec[:,0])),closphasdynspec[:,0] )
        plt.show()
        plt.scatter(np.arange(len(closphasdynspec[0,:])),closphasdynspec[0,:] )
        plt.show()
    
    fig = plt.figure(figsize=(5,5))
    plt.imshow(closphasdynspec, interpolation = 'none')
    plt.axis("off")
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.show()
    #plt.savefig('training_images/A/'+np.str(num)+'.png', box_inches='tight', dpi=64)
    plt.close()
    
    
    '''


rad2arcsec = 206265.
c = 3e8
doplot = 0
pix_size = 0.01
np.random.seed(1)
for i in range(4):
    randoms = np.random.rand(5,10000)
    N = 2084 # 0.1 arcsec each
    # put sources in central 20 pix (2 arcsec)
    x1 = (N/2-1/pix_size)+np.int(randoms[0,i]*(2/pix_size))
    y1 = (N/2-1/pix_size)+np.int(randoms[1,i]*(2/pix_size))
    x2 = (N/2-1/pix_size)+np.int(randoms[2,i]*(2/pix_size))
    y2 = (N/2-1/pix_size)+np.int(randoms[3,i]*(2/pix_size))
    dec = 60#randoms[4,i]*90



  
    make_sources(x1,y1,x2,y2,dec,i)
    
    


