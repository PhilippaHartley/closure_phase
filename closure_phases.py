from scipy.fftpack import ifftn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import sin, cos, sqrt, atan2, radians
from PIL import Image

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
    uv_centre = int(N/2)
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
    lamrang = np.linspace(1,1.05,100)
    hrang =  np.arange(0,120,0.4)

    dynspec = np.zeros((len(lamrang), len(hrang)))
    print (dynspec.shape)
    phase_plot = np.angle(Z)
    for lamcount, lam in enumerate(lamrang):
        
        # hour angle begins from same value each time, as if the sources all have same RA
        for hcount,h in enumerate(hrang):
            h = h*rad
            u = np.sin(h)*dx+np.cos(h)*dy
            v = -np.sin(raddec)*np.cos(h)*dx+np.sin(raddec)*np.sin(h)*dy+np.cos(raddec)*dz
        
       
            u/=(10000/100)
            v/=(10000/100)
            u*=20
            v*=20
            u = (u/np.float(lam)).astype(np.int)+uv_centre
            v = (v/np.float(lam)).astype(np.int)+uv_centre
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
    g = mkgauss([101,101], [50,50], 1, 50)


 

    blc1 = [int(x1-source_box_size[0]/2),int(y1-source_box_size[1]/2)]
    trc1 = [int(x1+source_box_size[0]/2),int(y1+source_box_size[1]/2)]
    blc2 = [int(x2-source_box_size[0]/2),int(y2-source_box_size[1]/2)]
    trc2 = [int(x2+source_box_size[0]/2),int(y2+source_box_size[1]/2)]



    xf[blc1[0]:trc1[0],blc1[1]:trc1[1]]+=g
    xf[blc2[0]:trc2[0],blc2[1]:trc2[1]]+=g
    #xf[x1, y1] = 10
    #xf[x1, y1] = 10

    plt.imshow(xf)
    plt.show()
    Z = np.fft.fftshift(np.fft.fft2(xf))
    plt.imshow(np.angle(Z))
    plt.show()

    plt.imshow((np.abs(Z)))
    plt.show()

    
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

    #a source at dec 90 will trace the most circular ellipse in uv



    # approximate radius of earth in km
    R = 6373.0
    
    lat1 = (lat_longs[2,0])*rad
    lon1 = (lat_longs[2,1])*rad
    lat2 = (lat_longs[0,0])*rad
    lon2 = (lat_longs[0,1])*rad
    lat3 = (lat_longs[7,0])*rad
    lon3 = (lat_longs[7,1])*rad
    
 

    p1,us1,vs1 = make_dynspec(lat1,lon1,lat2,lon2,Z)
    p2,us2,vs2 = make_dynspec(lat3,lon3,lat1,lon1,Z) 
    p3,us3,vs3 = make_dynspec(lat2,lon2,lat3,lon3,Z)


    plt.imshow(p1)
    plt.show()
    plt.imshow(p2)
    plt.show()
    plt.imshow(p3)
    plt.show()
    closphasdynspec = p1+p2-p3 

    print (us1, vs1)
    print (us2, vs2)
    print (us3, vs3)


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
    
    
    

    all_us = np.append(us1,us2)
    all_vs = np.append(vs1,vs2)
    all_us = np.append(all_us,us3)
    all_vs = np.append(all_vs,vs3)


    
    
np.random.seed(1)
for i in range(1):
    randoms = np.random.rand(5,10000)
    N = 513*2
    x1 = 20+np.int(randoms[0,i]*(N-40))
    y1 = 20+np.int(randoms[1,i]*(N-40))
    x2 = 20+np.int(randoms[2,i]*(N-40))
    y2 = 20+np.int(randoms[3,i]*(N-40))
    dec = 60#randoms[4,i]*90
    x1 = 80*2
    y1 = 250*2
    x2 = 350*2
    y2 = 250*2

  
    make_sources(x1,y1,x2,y2,dec,i)
    
    


