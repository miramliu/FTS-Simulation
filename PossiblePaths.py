from RayTraceFunctions import *
from BackgroundValues import *

'''Reflecting off of an ellipse'''
def ReflEll(Ray,thetL,originL,coeffellipse,center,ranges):
    if Ray is None:
        return
    Ray_Refl = []
    originG = [0,0,0] # the global origin
    thetG = [0,0,0] # rotation with respect to itself aka 0,0,0
    sourcepoint = Ray[2] #originalpoint
    v = Ray[3] #vector
    SPLi,VPLi = RT(sourcepoint,v,thetG,originG,thetL,originL) #point and vector in local coordinates
    pointsf,vectsf = REPCNi(coeffellipse,SPLi,VPLi) #IN LOCAL COORDINATES, REFLECTION OFF ELLIPSOID
    if SR3B(ranges, pointsf[0],pointsf[1],pointsf[2], center) == True:
        SPLf,VPLf = RT(pointsf,vectsf,thetL,originL,thetG,originG) #IN LOCAL COORDINATES, 
        Ray_Refl.append(Ray[0] + np.pi)
        Ray_Refl.append(Ray[1])
        Ray_Refl.append(SPLf)
        Ray_Refl.append(VPLf)
        Df = dist(sourcepoint,SPLf) #SO THE DISTANCE BETWEEN ORIGINAL POINT AND INTERSECTION OF ELLIPSOID.
        Ray_Refl.append(Ray[4] + Df) #should fuckin work.
        return Ray_Refl
    else:
        return

''' Entire run of the simulation (from input ray to the output ellipsoid, and all those below require the position of the mirror (as an origin). Only the eight paths that reach the detector are included, with each function referring to fork in the path chosen. For example, TTTTioM is the path of the ray that was transmitted through all four polarizers with the mirror at the position 'originM', and RTTRioM is the path of the ray that was reflected from polarizer 1, transmitted through polarizer 2 and 3, and and reflected from polarizer 4 with the mirror at the position 'originM'.'''

def TTTTioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)#OFF E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RRRRioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TTRRioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)    
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RTTRioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6) #E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RTRTioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5) #E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TRRTioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RRTTioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TRTRioM(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72


'''the functions below are identical to above, except are 'To Pickle' such that it returns every single ray created at every single step. All those below require the position of the mirror (as an origin)'''
def TTTTioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)#OFF E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri,Ray1,Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4,Ray_E72

def RRRRioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri,Ray1,Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4,Ray_E72

def TTRRioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)    
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4,Ray_E72

def RTTRioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6) #E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4, Ray_E72

def RTRTioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5) #E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4,Ray_E72

def TRRTioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4, Ray_E72

def RRTTioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originM) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4,Ray_E72

def TRTRioMPickle(Ri,p1,p2,p3,p4,originM):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originM) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ri, Ray1, Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4,Ray_E72




#all things commented out use RefEll_original (RefEllO)
#USING JUST TO DRAW OUT ALL THE PATHS POSSIBLE

'''
def TTTT(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4
    #return Ray_TP4
def RRRR(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4
def TTTR(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_RP4
def TTRT(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_TP4
def TRTT(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_TP4
def RTTT(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_TP4
def TTRR(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4
def TRTR(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4
def RTRT(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4
def TRRT(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4
def RRTT(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4
def RTTR(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4
def TRRR(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEllO(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_RP4
def RTRR(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEllO(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEllO(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_RP4
def RRTR(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEllO(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E5,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_RP4
def RRRT(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEllO(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEllO(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEllO(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEllO(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E6,coeffpolar,originpolar4,p4)
    return Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_TP4


#END OF DRAFTING/GRAPHING FUNCTIONS

def detectypos(Ri,p1,p2,p3,p4):
    R1 = TTTR(Ri,p1,p2,p3,p4)
    R2 = TTRT(Ri,p1,p2,p3,p4)
    R3 = TRTT(Ri,p1,p2,p3,p4)
    R4 = RTTT(Ri,p1,p2,p3,p4)
    R5 = TRRR(Ri,p1,p2,p3,p4)
    R6 = RTRR(Ri,p1,p2,p3,p4)
    R7 = RRTR(Ri,p1,p2,p3,p4)
    R8 = RRRT(Ri,p1,p2,p3,p4)
    return(R1,R2,R3,R4,R5,R6,R7,R8)

def detectyneg(Ri,p1,p2,p3,p4):
    R1 = TTTT(Ri,p1,p2,p3,p4)
    R2 = RRRR(Ri,p1,p2,p3,p4)
    R3 = TTRR(Ri,p1,p2,p3,p4)
    R4 = RTTR(Ri,p1,p2,p3,p4)
    R5 = RTRT(Ri,p1,p2,p3,p4)
    R6 = TRRT(Ri,p1,p2,p3,p4)
    R7 = RRTT(Ri,p1,p2,p3,p4)
    R8 = TRTR(Ri,p1,p2,p3,p4)
    Ray_E71 = ReflEll(R1,thet7,origin7,coeffellipse7)
    Ray_E72 = ReflEll(R2,thet7,origin7,coeffellipse7)
    Ray_E73 = ReflEll(R3,thet7,origin7,coeffellipse7)
    Ray_E74 = ReflEll(R4,thet7,origin7,coeffellipse7)
    Ray_E75 = ReflEll(R5,thet7,origin7,coeffellipse7)
    Ray_E76 = ReflEll(R6,thet7,origin7,coeffellipse7)
    Ray_E77 = ReflEll(R7,thet7,origin7,coeffellipse7)
    Ray_E78 = ReflEll(R8,thet7,origin7,coeffellipse7)
    return Ray_E71,Ray_E72,Ray_E73,Ray_E74,Ray_E75,Ray_E76,Ray_E77,Ray_E78

all the below are the entire fts run through given an initial ray and four polarizers 

def TTTTE(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def RRRRE(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def TTRRE(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)    
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def RTTRE(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def RTRTE(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def TRRTE(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def RRTTE(Ri,p1,p2,p3,p4):
    Ray_RP1 = IntPolR2(Ri,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7)
    return Ray_E72

def TRTRE(Ri,p1,p2,p3,p4):
    Ray_TP1 = IntPolT2(Ri,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7)
    return Ray_E72

all those below are ALL RAYS AT EVERY INTERSECTION

def TTTTioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)#OFF E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4,Ray_E72

def RRRRioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1, Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4,Ray_E72

def TTRRioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)    
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_TP1,Ray_E8,Ray_TP2,Ray_E3,Ray_M0,Ray_E4,Ray_RP3,Ray_E6,Ray_RP4,Ray_E72

def RTTRioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6) #E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4,Ray_E72

def RTRTioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5) #E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_RP1,Ray_E9,Ray_TP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4,Ray_E72

def TRRTioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_RP3,Ray_E5,Ray_TP4,Ray_E72

def RRTTioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_RP1,Ray_E9,Ray_RP2,Ray_E3,Ray_M0,Ray_E4,Ray_TP3,Ray_E5,Ray_TP4,Ray_E72

def TRTRioT(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray1,Ray_TP1,Ray_E8,Ray_RP2,Ray_E1,Ray_M0,Ray_E2,Ray_TP3,Ray_E6,Ray_RP4,Ray_E72


ALL those below are SOLELY input (from the source) and output.

def TTTTio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)#OFF E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RRRRio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TTRRio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_TP2 = IntPolT2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_TP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_RP3 = IntPolR2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_RP3, thet6,origin6,coeffellipse56,center6,range6) #off E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)    
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RTTRio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6) #E6
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RTRTio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E9
    Ray_TP2 = IntPolT2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_TP2,thet,origin1,coeffellipse,center1,range1) #E1
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E2
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5) #E5
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TRRTio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_RP3 = IntPolR2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_RP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def RRTTio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_RP1 = IntPolR2(Ray1,coeffpolar,originpolar1,p1)#p1
    Ray_E9 = ReflEll(Ray_RP1,thet5,origin9,coeffellipse56,center9,range9) #E8
    Ray_RP2 = IntPolR2(Ray_E9,coeffpolar,originpolar2,p2) #P2
    Ray_E3 = ReflEll(Ray_RP2,thet,origin3,coeffellipse,center3,range3) #E3
    Ray_M0 = IntM2(Ray_E3, coeffmirr, originG) #off mirror
    Ray_E4 = ReflEll(Ray_M0, thet,origin4,coeffellipse,center4,range4) #off E4
    Ray_TP3 = IntPolT2(Ray_E4,coeffpolar,originpolar3,p3) #P3
    Ray_E5 = ReflEll(Ray_TP3, thet5,origin5,coeffellipse56,center5,range5)
    Ray_TP4 = IntPolT2(Ray_E5,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_TP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72

def TRTRio(Ri,p1,p2,p3,p4):
    Ray1 = ReflEll(Ri,thet10,origin10,coeffellipse7,center10,range10)
    Ray_TP1 = IntPolT2(Ray1,coeffpolar,originpolar1,p1) #P1
    Ray_E8 = ReflEll(Ray_TP1,thet6,origin8,coeffellipse56,center8,range8) #E8
    Ray_RP2 = IntPolR2(Ray_E8,coeffpolar,originpolar2,p2) #P2
    Ray_E1 = ReflEll(Ray_RP2,thet,origin1,coeffellipse,center1,range1) #E3
    Ray_M0 = IntM2(Ray_E1, coeffmirr, originG) #off mirror
    Ray_E2 = ReflEll(Ray_M0, thet,origin2,coeffellipse,center2,range2) #off E4
    Ray_TP3 = IntPolT2(Ray_E2,coeffpolar,originpolar3,p3) #P3
    Ray_E6 = ReflEll(Ray_TP3, thet6,origin6,coeffellipse56,center6,range6)
    Ray_RP4 = IntPolR2(Ray_E6,coeffpolar,originpolar4,p4)
    Ray_E72 = ReflEll(Ray_RP4,thet7,origin7,coeffellipse7,center7,range7)
    return Ray_E72'''