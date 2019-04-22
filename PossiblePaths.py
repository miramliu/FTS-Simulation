'''These are the explicit functions of the specific paths of rays (transmission and reflection through four polarizers and off of all the mirrors/ellipsoids) through the simulation of the fourier transform spectrometer. Mira Liu 04/21/19  '''

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
