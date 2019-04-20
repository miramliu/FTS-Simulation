README

This is a repository of a version of the ray trace simulation of a compact fourier transform spectrometer written by Mira Liu over the years 2016-2019 with Dr. Stephan S. Meyer, Dr. Ritoban Basu Thakur, and Zhaodi Pan. 






Output of files ending in "_topickle"
- these were written to explicitly save every aspect of the ray trace simulation by following and saving each ray. 
# summary of the array "Rays"
### RAYS ARE ORGANIZED AS FOLLOWS

RAYS[1stRay
     [Mirror Position1
      [[Path1
       [[section1],[section2],...[section12]]
       ]
      [Path2
       [[section1],...[section12]]
      ]
       ...
      [Path8
          
      ]
     ]
     Mirror Position2
      [[path1[sections]]
       [path2[sections]]
       ...
       [path8[sections]]
      ]
     ]   
2ndRay[Mirror Positions[Paths[sections]]]
...
NthRay[Mirror Positions[Paths[sections]]]]

So generally Rays[Ns[Mirror Positions[Paths[sections]]]]

Mirror position is -18 to 18 divided into Nsize segments (Nsize being number of samples taken as the mirror travels)

so the Mirror Position = "7" is the location of the mirror at the 7th position of the mirror travelling from -18 to 18

path is which of the 8 different combinations of transmission and reflection that lead to the detector

in order 0-7:TTTT,RRRR,TTRR,RTTR,RTRT,TRRT,RRTT,TRTR, with T=transmission, R=reflection

there are 12 sections of each path (from launch from source to launch from last ellipsoid)

quick tutorial:
if you want a certain ray "n"
Rays[n] will return the Ray Number, ALL mirror positions, paths, and sections of that one ray.
Rays[n][0] will return the name "Ray: n"
Rays[n][1] will return the first mirror position, corresponding paths, and sections.
Rays[n][2] will return the second mirror positions,corresponding paths,and sections
Rays[n][m][0] will return the mth mirror position's actual location (x,y,z)
Rays[n][m][1] will return the first path and sections (TTTT and 12 corresponding sections)
Rays[n][m][2] will return the second path and sections (RRRR and 12 corresponding sections)
and of course each section is defined by a singly Ray of the form 
Ray = [polarization(theta), intensity (magnitude of vector), point (x,y,z), direction(Vx,Vy,Vz), distance travelled (D)]

So to get ALL the LAST sections of all Rays (all of the rays reflecting off of the last ellipsoid into the detector plane: 
for i in range(N): #number of rays
    for j in range(1,Nsize): #mirror position (0 is the name of the mirror position)
        for k in range(1,9): #which paths (0 is the name of the path)
            print(i,j,k,Rays[i][j][k][12]) #12th section



