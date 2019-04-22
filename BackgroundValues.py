from RayTraceFunctions import *

## the four polarizers in order.
p1 = np.pi/4
p2 = np.pi/2
p3 = np.pi
p4 = np.pi/4

originG = [0.,0.,0.] # the global origin
thetG = [0.,0.,0.] # rotation with respect to itself aka 0,0,0
origin1 = [-32.075,-128.,0.] #x,y (ellipse1)
origin2 = [32.075,-128.,0.] #x,y (ellipse2)
origin3 = [-32.075,128.,0.] #x,y (ellipse3)
origin4 = [32.075,128,0.] #x,y  (ellipse4)
origin5 = [96.225,-120.50,0.] # (ellipse5)
origin6 = [96.225,120.50,0.] # (ellipse6)
origin7 = [128.3,7.5,39.850000] # (ellipse7)
origin8 = [-96.225,-120.50,0.] # (ellipse8)
origin9 = [-96.225,120.50,0.] # (ellipse9)
origin10 = [-128.3,7.5,-39.850000] # last ellipsoid

#coefficients are like [a,b,c] of an ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2
coeffellipse7 = [164.54585247700001,99.690818975602866,130.9086635] #for ellipse 7
coeffellipse56 = [256.65344272795591,248.39387505453516,64.58693753]  #for ellipses 5&6&8&9
coeffellipse = [263.915180503,256.0,64.15] #for the four center ellipses

coeffmirr = [31.75,25.4,19.05] #for the mirror
coeffpolar = [32.075,32.075,0] #for polarizers (2d circle)

thet = [0,0,0] 
thet5 = [0,0,-.116385] #angle of rotation
thet6 = [0,0,.116385] #angle of rotation
thet7 = [0,0.309319724356,1.31064594453] #angle of rotation
thet10 = [0,0.309319724356,-1.31064594453] #angle of rotation

originpolar1 = [-128.3,0.0,0.0] #global origin of polarizer 1
originpolar2 = [-64.15,0.0,0.0] #global origin of polarizer 2
originpolar3 = [64.15,0.0,0.0] #global origin of polarizer 3
originpolar4 = [128.3,0.0,0.0] #global origin of polarizer 4

origindet = [160.375,-113] #global origin of detector
radet = 7.9375 #radius of detector 

#below are the geometric centers and ranges of aspects of the FTS.
center1,range1= [0.0, 254.99883040161689, 0.0], [31.427020202020202, 200, 31.400097443962792]
center2,range2= [-3.5527136788005009e-15, 254.99883040161689, 0.0], [31.427020202020202, 200, 31.400097443962792]
center3,range3= [0.0, -254.99883040161689, 0.0], [31.427020202020202, 200, 31.400097443962792]
center4,range4= [-3.5527136788005009e-15, -254.99883040161689, 0.0], [31.427020202020202, 200, 31.400097443962792]
center5,range5= [28.333333333333336, 244.95832095066169, 0.0], [30.909090909090921, 200, 31.315767057290568]
center6,range6= [28.333333333333336, -244.95832095066169, 0.0], [30.909090909090921, 200, 31.315767057290568]
center8,range8= [-28.333333333333329, 244.95832095066169, 0.0], [30.909090909090914, 200, 31.315767057290568]
center9,range9= [-28.333333333333329, -244.95832095066169, 0.0], [30.909090909090914, 200, 31.315767057290568]
center10,range10= [-96.458686868686868, 3.279771448324329, 73.24759051407338], [41.577020202020229, 33.648217777841843, 20.001702052348694]
center7,range7= [96.45868686868684, 3.1907688264097978, -73.870253223692842], [41.577020202020208, 200, 19.642335814884394]
Ecenter7 = [192.45-32.075,0,0]
Ecenter10 = [-128.3-32.075,0,0]
