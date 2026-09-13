# Chamber A dictionary generator, derived from createSMTCdict.py (Mathias Gruber, 2012).
# Run in a new case directory. See ../README.md for usage.
# Import math package
from math import sin, cos, radians, asin, atan, sqrt, tan
from pathlib import Path
Path("system").mkdir(parents=True,exist_ok=True)
import re

###########################################
# Description                             #
# This script creates the blockMeshDict   #
# for creating the SMTC chamber using     #
# blockMesh. Various parameters can be    #
# adjusted using the script               #
# All units are in mm                     #
###########################################

###########################################
#       USER INPUT - START                #
###########################################

# Main Chamber
main_chamber_height = 1.;
main_chamber_length = 30.;
main_chamber_width = 15;
main_chamber_dist_to_inner = 1;
main_chamber_cornorRadii = 3.

# Inlet
inlet_cylinder_angle = 48.7;
inlet_cylinder_radius = 0.75;
inlet_cylinder_length = 10.;
inlet_cylinder_squareRadius = 0.25;
inlet_cylinder_MC_sep = 0.05;
inlet_cylinder_edgeDist = 1.6;

# Carvings - lower BOX
carvings_height = 1
carvings_length = 4.85
carvings_width = 3.
carvings_edgeDist = 1.5;
carvings_inner_center = 0.5;

# Refinement Option, only integers over or equal to 1.
mesh_Refinement = 1;

###########################################
#       Script - HEADER                   #
###########################################
# Open File
filename = "system/blockMeshDict";
FILE = open(filename,"w");

header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolboy           |
|  \\    /   O peration     | Version:  1.7.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// Dimensions in milli-meters
scale                   0.001;
""";

FILE.write(header);

MC_l = main_chamber_length/2
MC_h = main_chamber_height/2
MC_w = main_chamber_width/2
MC_sep = main_chamber_dist_to_inner;
MC_R = main_chamber_cornorRadii;
edge_factor = 0.75;

I_r = inlet_cylinder_radius;
I_a = inlet_cylinder_angle
I_cr = inlet_cylinder_squareRadius
I_l = inlet_cylinder_length
I_MC_sep = inlet_cylinder_MC_sep
I_e = inlet_cylinder_edgeDist

iter_moveCyl = 4

c_h = carvings_height
c_l = carvings_length
c_w = carvings_width
c_eD = carvings_edgeDist
c_iS = carvings_inner_center

atc_factor = 0;
atc1 = 0.060*atc_factor;
atc2 = 0.02*atc_factor;
atc3 = 0.019*atc_factor;
atc4 = 0.15*atc_factor;
atc5 = 0.4*atc_factor;
atc6 = 0.60*atc_factor;

MR = mesh_Refinement
MA = 1.1; # Adjusts the inlet/carving distance to remove artefacts

###########################################
#       Script - FUNCTIONS                #
###########################################
def sinD( x ):
    return sin( radians(x) );

def cosD( x ):
    return cos( radians(x) );

def tanD( x ):
    return tan( radians(x) );

def rot3d( oriX , oriY , oriZ , posX , posY , posZ , theta , phi , psi ):
    M11 = cosD(theta)*cosD(psi);
    M12 = -cosD(phi)*sinD(psi)+sinD(phi)*sinD(theta)*cosD(psi);
    M13 = sinD(phi)*sinD(psi) + cosD(phi)*sinD(theta)*cosD(psi);
    M21 = cosD(theta)*sinD(psi);
    M22 = cosD(phi)*cosD(psi)+sinD(phi)*sinD(theta)*sinD(psi);
    M23 = -sinD(phi)*cosD(psi)+cosD(phi)*sinD(theta)*sinD(psi);
    M31 = -sinD(theta);
    M32 = sinD(phi)*cosD(theta);
    M33 = cosD(phi)*cosD(theta);

    xNew = oriX + (posX-oriX)*M11 + (posY-oriY)*M12 + (posZ-oriZ)*M13;
    yNew = oriY + (posX-oriX)*M21 + (posY-oriY)*M22 + (posZ-oriZ)*M23;
    zNew = oriZ + (posX-oriX)*M31 + (posY-oriY)*M32 + (posZ-oriZ)*M33;

    return [xNew , yNew , zNew];

# Attach central inlet to main chamber, return [ y , z ], p.32 in notebook
def attach( posZ , posY , movement):
    # New Positions
    vector = rot3d( 0 , 0 , 0 , 0 , movement , 0 , 0 ,-I_a , 0 );
    print(vector);
    yNew = posY - vector[1]
    zNew = posZ + vector[2]

    # New xz positions
    return [ yNew , zNew ];

def attach2( posX ,  posY , posZ):
    # Determine move vector
    vector = rot3d( 0 , 0 , 0 , 0 , 0 , inlet_cylinder_MC_sep , 0 , -I_a , 0 );
    print(vector);
    # Determine distance
    distance = inlet_cylinder_MC_sep*2;
    eps = inlet_cylinder_MC_sep*MA; i = 0;
    while ( eps < distance and i < 100):
        posY = posY - vector[1]
        y_Carv = MC_l - c_eD
        distance = posY - y_Carv

        i = i +1 ;
        print(distance);

    # New xz positions
    return [ posX, posY , posZ ];

def attach3( posX ,  posY , posZ):
    # Determine move vector
    vector = rot3d( 0 , 0 , 0 , 0 , 0 , inlet_cylinder_MC_sep , 0 , -I_a , 0 );
    print(vector);

    # Determine distance
    distance = inlet_cylinder_MC_sep*2;
    eps = inlet_cylinder_MC_sep*MA; i = 0;
    while ( eps < distance and i < 100):
        posX = posX - vector[0]
        posY = posY - vector[1]
        posZ = posZ - vector[2]
        # Determine distance
        alpha = atan( (posZ-(MC_h+c_h)) / (posY-(MC_l-c_eD-c_w/2.)) );
        print("Angle found to be",alpha)
        circleZ_Y = [ sin(alpha)*(c_w/2.)+MC_h+c_h , MC_l - c_eD - c_w/2. + cos(alpha)*(c_w/2.) ]
        distance = sqrt((posZ - circleZ_Y[0])**2 + (posY - circleZ_Y[1])**2)

        i = i +1 ;
        print(distance);

    # New xz positions
    return [ posX, posY , posZ ];

###########################################
#       Script - Vertices                 #
###########################################
main_chamber_vertices = [];
main_chamber_vertices.append( [ MC_w                                    , 0                                         , 0 ]               ); # 0
main_chamber_vertices.append( [ MC_w                                    , 0                                         , MC_h ]            ); # 1
main_chamber_vertices.append( [ MC_w-MC_sep                             , 0                                         , 0 ]               ); # 2
main_chamber_vertices.append( [ MC_w-MC_sep                             , 0                                         , MC_h ]            ); # 3
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , 0                                         , 0 ]               ); # 4
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , 0                                         , MC_h ]            ); # 5
main_chamber_vertices.append( [ MC_w                                    , MC_l-2*MC_sep                             , 0 ]               ); # 6
main_chamber_vertices.append( [ MC_w                                    , MC_l-2*MC_sep                             , MC_h ]            ); # 7
main_chamber_vertices.append( [ MC_w-MC_sep                             , MC_l-2*MC_sep                             , 0 ]               ); # 8
main_chamber_vertices.append( [ MC_w-MC_sep                             , MC_l-2*MC_sep                             , MC_h ]            ); # 9
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l-2*MC_sep                             , 0 ]               ); # 10
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l-2*MC_sep                             , MC_h ]            ); # 11
main_chamber_vertices.append( [ MC_w-MC_sep+sinD(45)*MC_sep*edge_factor , MC_l-MC_sep+sinD(45)*MC_sep*edge_factor   , 0 ]               ); # 12
main_chamber_vertices.append( [ MC_w-MC_sep+sinD(45)*MC_sep*edge_factor , MC_l-MC_sep+sinD(45)*MC_sep*edge_factor   , MC_h ]            ); # 13
main_chamber_vertices.append( [ MC_w-MC_sep                             , MC_l-MC_sep                               , 0 ]               ); # 14
main_chamber_vertices.append( [ MC_w-MC_sep                             , MC_l-MC_sep                               , MC_h ]            ); # 15
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l-MC_sep                               , 0 ]               ); # 16
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l-MC_sep                               , MC_h ]            ); # 17
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l                                      , 0 ]               ); # 18
main_chamber_vertices.append( [ MC_w-MC_sep*2                           , MC_l                                      , MC_h ]            ); # 19
main_chamber_vertices.append( [ 0                                       , 0                                         , 0 ]               ); # 20
main_chamber_vertices.append( [ 0                                       , 0                                         , MC_h ]            ); # 21
main_chamber_vertices.append( [ 0                                       , MC_l-2*MC_sep                             , 0 ]               ); # 22
main_chamber_vertices.append( [ 0                                       , MC_l-2*MC_sep                             , MC_h ]            ); # 23
main_chamber_vertices.append( [ 0                                       , MC_l-MC_sep                               , 0 ]               ); # 24
main_chamber_vertices.append( [ 0                                       , MC_l-MC_sep                               , MC_h ]            ); # 25
main_chamber_vertices.append( [ 0                                       , MC_l                                      , 0 ]               ); # 26
main_chamber_vertices.append( [ 0                                       , MC_l                                      , MC_h ]            ); # 27

inlet_vertices = [];
inlet_vertices.append(        [ 0                                       , MC_l - 2*I_r - I_e                             , MC_h ]            ); # 28 0
inlet_vertices.append(        [ 0                                       , MC_l - 2*I_r - I_e                              , MC_h + I_l ]      ); # 29 1
inlet_vertices.append(        [ 0                                       , MC_l - I_r - I_cr - I_e                         , MC_h ]            ); # 30 2
inlet_vertices.append(        [ 0                                       , MC_l - I_r - I_cr - I_e                         , MC_h + I_l ]      ); # 31 3
inlet_vertices.append(        [ 0                                       , MC_l - I_r + I_cr - I_e                         , MC_h ]            ); # 32 4
inlet_vertices.append(        [ 0                                       , MC_l - I_r + I_cr - I_e                         , MC_h + I_l ]      ); # 33 5
inlet_vertices.append(        [ 0                                       , MC_l - I_e                                      , MC_h ]            ); # 34 6
inlet_vertices.append(        [ 0                                       , MC_l - I_e                                      , MC_h + I_l ]      ); # 35 7
inlet_vertices.append(        [ I_cr                                    , MC_l - I_r - I_cr - I_e                         , MC_h ]            ); # 36 8
inlet_vertices.append(        [ I_cr                                    , MC_l - I_r - I_cr - I_e                         , MC_h + I_l ]      ); # 37 9
inlet_vertices.append(        [ I_cr                                    , MC_l - I_r + I_cr - I_e                         , MC_h ]            ); # 38 10
inlet_vertices.append(        [ I_cr                                    , MC_l - I_r + I_cr - I_e                         , MC_h + I_l ]      ); # 39 11
inlet_vertices.append(        [ I_r*sinD(45)                            , MC_l - I_r - sinD(45)*I_r - I_e                 , MC_h ]            ); # 40 12
inlet_vertices.append(        [ I_r*sinD(45)                            , MC_l - I_r - sinD(45)*I_r - I_e                 , MC_h + I_l ]      ); # 41 13
inlet_vertices.append(        [ I_r*sinD(45)                            , MC_l - I_r + sinD(45)*I_r - I_e                 , MC_h ]            ); # 42 14
inlet_vertices.append(        [ I_r*sinD(45)                            , MC_l - I_r + sinD(45)*I_r - I_e                 , MC_h + I_l ]      ); # 43 15

carving_vertices = [];
# Carving bottom box, page 42
carving_vertices.append(      [ c_l                                     , MC_l - c_w - c_eD                                         , MC_h ]   ); # 44
carving_vertices.append(      [ c_l                                     , MC_l - c_w - c_eD                                         , MC_h+c_h   ]      ); # 45
carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h ]   ); # 46
carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 47
carving_vertices.append(      [ 0                                       , MC_l - c_w - c_eD                                         , MC_h ]   ); # 48
carving_vertices.append(      [ 0                                       , MC_l - c_w - c_eD                                         , MC_h+c_h   ]      ); # 49
carving_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h ]   ); # 50
carving_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 51
carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h ]   ); # 52
carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 53
carving_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h ]   ); # 54
carving_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 55
carving_vertices.append(      [ c_l                                     , MC_l              - c_eD                                  , MC_h ]   ); # 56
carving_vertices.append(      [ c_l                                     , MC_l              - c_eD                                  , MC_h+c_h   ]      ); # 57
carving_vertices.append(      [ 0                                       , MC_l              - c_eD                                  , MC_h ]   ); # 58
carving_vertices.append(      [ 0                                       , MC_l              - c_eD                                  , MC_h+c_h   ]      ); # 59
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h ]   ); # 60
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 61
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h ]   ); # 62
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 63
carving_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD - cosD(45)*(c_w*0.5)             , MC_h ]   ); # 64
carving_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD - cosD(45)*(c_w*0.5)             , MC_h+c_h   ]      ); # 65
carving_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD + cosD(45)*(c_w*0.5)             , MC_h ]   ); # 66
carving_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD + cosD(45)*(c_w*0.5)             , MC_h+c_h   ]      ); # 67

carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 68
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 69
carving_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 70
carving_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 71

carving_vertices.append( rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , 0) ); # 72
carving_vertices.append( rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 135 , 0) ); # 73
carving_vertices.append( rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , -45) ); # 74
carving_vertices.append( rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , -135) ); # 75

carving_vertices.append(      [ 0                                     , MC_l              - c_eD                                  , MC_h+c_h   ]      ); # 76
carving_vertices.append(      [ 0                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 77
carving_vertices.append(      [ 0                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h   ]      ); # 78
carving_vertices.append( rot3d( 0     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   0     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , 0) ); # 79
carving_vertices.append(      [ 0                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 80
carving_vertices.append(      [ 0                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , MC_h+c_h+c_iS   ]      ); # 81
carving_vertices.append( rot3d( 0     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   0     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 135 , 0) ); # 82
carving_vertices.append(      [ 0                                     , MC_l - c_w - c_eD                                         , MC_h+c_h   ]      ); # 83

# New main chamber vertices
Newmain_chamber_vertices = [];


Newmain_chamber_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , 0 ]   ); # 84
Newmain_chamber_vertices.append(      [ 0                                       , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , 0 ]   ); # 85
Newmain_chamber_vertices.append(      [ 0                                       , MC_l              - c_eD                                  , 0 ]   ); # 86
Newmain_chamber_vertices.append(      [ c_l                                     , MC_l - c_w - c_eD                                         , 0 ]   ); # 87
Newmain_chamber_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , 0 ]   ); # 88
Newmain_chamber_vertices.append(      [ c_l                                     , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , 0 ]   ); # 89
Newmain_chamber_vertices.append(      [ c_l                                     , MC_l              - c_eD                                  , 0 ]   ); # 90
Newmain_chamber_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)-c_iS/2 - c_eD                           , 0 ]   ); # 91
Newmain_chamber_vertices.append(      [ c_l+c_iS                                , MC_l - c_w*(1/2.)+c_iS/2 - c_eD                           , 0 ]   ); # 92
Newmain_chamber_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD - cosD(45)*(c_w*0.5)             , 0 ]   ); # 93
Newmain_chamber_vertices.append(      [ c_l+sinD(45)*(c_w*0.5)                  , MC_l - c_w*(1./2) - c_eD + cosD(45)*(c_w*0.5)             , 0 ]   ); # 94

Newmain_chamber_vertices.append(      [ MC_w                                     , MC_l - MC_R*1.5                                               , 0 ]   ); # 95
Newmain_chamber_vertices.append(      [ MC_w                                     , MC_l - MC_R*1.5                                               , MC_h ]   ); # 96
Newmain_chamber_vertices.append(      [ MC_w-MC_R+MC_R*sinD(45)                  , MC_l-MC_R+MC_R*sinD(45)                                   , 0 ]   ); # 97
Newmain_chamber_vertices.append(      [ MC_w-MC_R+MC_R*sinD(45)                  , MC_l-MC_R+MC_R*sinD(45)                                   , MC_h ]   ); # 98
Newmain_chamber_vertices.append(      [ MC_w-MC_R                                , MC_l                                                      , 0 ]   ); # 99
Newmain_chamber_vertices.append(      [ MC_w-MC_R                                , MC_l                                                      , MC_h ]   ); # 100
Newmain_chamber_vertices.append(      [ 0                                       , MC_l - c_w - c_eD                                         , 0 ]   ); # 101

# Vertice List
vertices = [];

# Main Chamber
i = 0;
for vertex in main_chamber_vertices:
    vertices.append( "    ( "+str(vertex[0])+"    "+str(vertex[1])+"    "+str(vertex[2])+"    )   // Number "+str(i) );
    i = i + 1;

# Inlet Chamber
i = 0;
for point in inlet_vertices:
    inlet_vertices[i] = rot3d( 0 , MC_l , MC_h , inlet_vertices[i][0],  inlet_vertices[i][1] ,  inlet_vertices[i][2] , 0 , -I_a , 0);
    i = i + 1;

[ inlet_vertices[0][0] , inlet_vertices[0][1] , inlet_vertices[0][2] ]   = attach3( inlet_vertices[0][0] ,inlet_vertices[0][1] ,inlet_vertices[0][2] )
[ inlet_vertices[2][0] , inlet_vertices[2][1] , inlet_vertices[2][2] ]   = attach3( inlet_vertices[2][0] ,inlet_vertices[2][1] ,inlet_vertices[2][2] )
[ inlet_vertices[8][0] , inlet_vertices[8][1] , inlet_vertices[8][2] ]   = attach3( inlet_vertices[8][0] ,inlet_vertices[8][1] ,inlet_vertices[8][2] )
[ inlet_vertices[4][0] , inlet_vertices[4][1] , inlet_vertices[4][2] ]   = attach3( inlet_vertices[4][0] ,inlet_vertices[4][1] ,inlet_vertices[4][2] )
[ inlet_vertices[10][0] , inlet_vertices[10][1] , inlet_vertices[10][2] ]   = attach3( inlet_vertices[10][0] ,inlet_vertices[10][1] ,inlet_vertices[10][2] )

[ inlet_vertices[12][0] , inlet_vertices[12][1] , inlet_vertices[12][2] ]   = attach3( inlet_vertices[12][0] ,inlet_vertices[12][1] ,inlet_vertices[12][2] )
[ inlet_vertices[14][0] , inlet_vertices[14][1] , inlet_vertices[14][2] ]     = attach2( inlet_vertices[14][0], inlet_vertices[14][1], inlet_vertices[14][2] )
[ inlet_vertices[6][0] , inlet_vertices[6][1] , inlet_vertices[6][2] ]     = attach2( inlet_vertices[6][0], inlet_vertices[6][1], inlet_vertices[6][2] )

i = 0;
for point in inlet_vertices:
    if inlet_vertices[i][2] < I_l*sinD(I_a):
        a = 2;#[ inlet_vertices[i][1] , inlet_vertices[i][2] ]     = attach( inlet_vertices[i][2], inlet_vertices[i][1] )
    i = i + 1;

i = 28;
for vertex in inlet_vertices:
    vertices.append( "    ( "+str(vertex[0])+"    "+str(vertex[1])+"    "+str(vertex[2])+"    )   // Number "+str(i) );
    i = i + 1;

# Carving


i = 44;
for vertex in carving_vertices:
    vertices.append( "    ( "+str(vertex[0])+"    "+str(vertex[1])+"    "+str(vertex[2])+"    )   // Number "+str(i) );
    i = i + 1;

# New main chamber
i = 87;
for vertex in Newmain_chamber_vertices:
    vertices.append( "    ( "+str(vertex[0])+"    "+str(vertex[1])+"    "+str(vertex[2])+"    )   // Number "+str(i) );
    i = i + 1;


# Add points to buffer
buffer = """
// Points in {xyz}
vertices
(
""";
for vertice in vertices:
    buffer = buffer+ vertice+"\n";


# Close buffer
buffer = buffer+");\n";

# Write to file
FILE.write(buffer);


###########################################
#       Script - Blocks                    #
###########################################
buffer = "";
buffer = """
// Blocks
blocks
(
""";

blocks = [];

# Inlet
blocks.append( "    hex ( 28 40 36 30 29 41 37 31 )       ( "+str(MR*3)+" "+str(MR*4)+" "+str(MR*20)+" ) simpleGrading ( 1 1 5 )"   ); # Blok 9
blocks.append( "    hex ( 36 40 42 38 37 41 43 39 )       ( "+str(MR*4)+" "+str(MR*5)+" "+str(MR*20)+" ) simpleGrading ( 1 1 5 )"   ); # Blok 10
blocks.append( "    hex ( 32 38 42 34 33 39 43 35 )       ( "+str(MR*3)+" "+str(MR*4)+" "+str(MR*20)+" ) simpleGrading ( 1 1 5 )"   ); # Blok 11
blocks.append( "    hex ( 30 36 38 32 31 37 39 33 )       ( "+str(MR*3)+" "+str(MR*5)+" "+str(MR*20)+" ) simpleGrading ( 1 1 5 )"   ); # Blok 12

# Carving
blocks.append( "    hex ( 48 44 46 50 49 45 47 51 )       ( "+str(MR*30)+" "+str(MR*5)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 13
blocks.append( "    hex ( 50 46 52 54 51 47 53 55 )       ( "+str(MR*30)+" "+str(MR*10)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 14
blocks.append( "    hex ( 54 52 56 58 55 53 57 59 )       ( "+str(MR*30)+" "+str(MR*5)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 15
blocks.append( "    hex ( 44 64 60 46 45 65 61 47 )       ( "+str(MR*5)+" "+str(MR*5)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 16
blocks.append( "    hex ( 60 64 66 62 61 65 67 63 )       ( "+str(MR*5)+" "+str(MR*10)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 17
blocks.append( "    hex ( 52 62 66 56 53 63 67 57 )       ( "+str(MR*5)+" "+str(MR*5)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 18
blocks.append( "    hex ( 46 60 62 52 47 61 63 53 )       ( "+str(MR*5)+" "+str(MR*10)+" "+str(MR*10)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 19

blocks.append( "    hex ( 53 47 61 63 71 68 69 70 )       ( "+str(MR*10)+" "+str(MR*5)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 20
blocks.append( "    hex ( 57 53 63 67 72 71 70 74 )       ( "+str(MR*5)+" "+str(MR*5)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 21
blocks.append( "    hex ( 67 63 61 65 74 70 69 75 )       ( "+str(MR*5)+" "+str(MR*10)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 22
blocks.append( "    hex ( 47 45 65 61 68 73 75 69 )       ( "+str(MR*5)+" "+str(MR*5)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 23
blocks.append( "    hex ( 71 68 69 70 72 73 75 74 )       ( "+str(MR*10)+" "+str(MR*5)+" "+str(MR*5)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 24
blocks.append( "    hex ( 53 55 51 47 71 80 81 68 )       ( "+str(MR*30)+" "+str(MR*10)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 25
blocks.append( "    hex ( 47 51 49 45 68 81 82 73  )       ( "+str(MR*30)+" "+str(MR*5)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 26
blocks.append( "    hex ( 71 80 81 68 72 79 82 73  )       ( "+str(MR*30)+" "+str(MR*10)+" "+str(MR*5)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 27
blocks.append( "    hex ( 57 59 55 53 72 79 80 71  )       ( "+str(MR*30)+" "+str(MR*5)+" "+str(MR*8)+" ) simpleGrading ( 1 1 1 )"   ); # Blok 28

# Main Chamber Blocks - NEW
blocks.append( "    hex ( 101 87 88 84 48 44 46 50 )       ( "+str(MR*30)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 29
blocks.append( "    hex ( 84 88 89 85 50 46 52 54 )        ( "+str(MR*30)+" "+str(MR*10)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 30
blocks.append( "    hex ( 85 89 90 86 54 52 56 58 )        ( "+str(MR*30)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 31
blocks.append( "    hex ( 87 93 91 88 44 64 60 46 )        ( "+str(MR*5)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 32
blocks.append( "    hex ( 88 91 92 89 46 60 62 52 )        ( "+str(MR*5)+" "+str(MR*10)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 33
blocks.append( "    hex ( 89 92 94 90 52 62 66 56 )        ( "+str(MR*5)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 34
blocks.append( "    hex ( 91 93 94 92 60 64 66 62 )        ( "+str(MR*5)+" "+str(MR*10)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 35
blocks.append( "    hex ( 86 90 99 26 58 56 100 27 )        ( "+str(MR*30)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 36
blocks.append( "    hex ( 90 94 97 99 56 66 98 100 )        ( "+str(MR*5)+" "+str(MR*5)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 37
blocks.append( "    hex ( 93 95 97 94 64 96 98 66 )        ( "+str(MR*5)+" "+str(MR*10)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 38
blocks.append( "    hex ( 2 0 95 93 3 1 96 64 )        ( "+str(MR*5)+" "+str(MR*30)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 39
blocks.append( "    hex ( 4 2 93 87 5 3 64 44 )        ( "+str(MR*5)+" "+str(MR*30)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 40
blocks.append( "    hex ( 20 4 87 101 21 5 44 48 )        ( "+str(MR*30)+" "+str(MR*30)+" 30 ) simpleGrading ( 1 1 5 )"   ); # Blok 41

# Add points to buffer
for block in blocks:
    buffer = buffer+ block+"\n";

# Close buffer
buffer = buffer+");\n";

# Write to file
FILE.write(buffer);

###########################################
#       Script - Edges                    #
###########################################

print("\nCalculating edge arc positions ");

buffer = "";
buffer = """
// Edges
edges
(
""";
edges = [];

### MAIN CHAMBER ###
MC_edges = []
MC_edges.append( [ "6 12"  , MC_w - MC_sep + cosD(45./2)*MC_sep*edge_factor           , MC_l -MC_sep + sinD(45./2)*MC_sep*edge_factor , 0           ] )
MC_edges.append( [ "7 13"  , MC_w - MC_sep + cosD(45./2)*MC_sep*edge_factor           , MC_l -MC_sep + sinD(45./2)*MC_sep*edge_factor , MC_h        ] )
MC_edges.append( [ "12 18" , MC_w - MC_sep + sinD(45./2)*MC_sep*edge_factor           , MC_l -MC_sep + cosD(45./2)*MC_sep*edge_factor , 0           ] )
MC_edges.append( [ "13 19" , MC_w - MC_sep + sinD(45./2)*MC_sep*edge_factor           , MC_l -MC_sep + cosD(45./2)*MC_sep*edge_factor , MC_h        ] )

I_edges = []
I_edges.append( [ "28 40" , sinD(45/2)*I_r                                           , MC_l - I_r - cosD(45/2)*I_r - I_e                   , MC_h        ] ) # 0
I_edges.append( [ "29 41" , sinD(45/2)*I_r                                           , MC_l - I_r - cosD(45/2)*I_r - I_e                   , MC_h + I_l  ] ) # 1
I_edges.append( [ "40 42" , I_r                                                      , MC_l - I_r - I_e                                    , MC_h        ] ) # 2
I_edges.append( [ "41 43" , I_r                                                      , MC_l - I_r - I_e                                    , MC_h + I_l  ] ) # 3
I_edges.append( [ "42 34" , sinD(45/2)*I_r                                           , MC_l - I_r + cosD(45/2)*I_r - I_e                   , MC_h        ] ) # 4
I_edges.append( [ "43 35" , sinD(45/2)*I_r                                           , MC_l - I_r + cosD(45/2)*I_r - I_e                   , MC_h + I_l  ] ) # 5

## MAIN CHAMBER
i = 0;
for edge in MC_edges:
    edges.append( "  arc  "+MC_edges[i][0]+"    ( "+str(MC_edges[i][1])+"  "+str(MC_edges[i][2])+"  "+str(MC_edges[i][3])+" )" );
    i = i + 1;

## INLET
i = 0;
for edge in I_edges:
    [ I_edges[i][1] , I_edges[i][2] , I_edges[i][3] ] = rot3d( 0 , MC_l , MC_h , I_edges[i][1],  I_edges[i][2] ,  I_edges[i][3] , 0 , -I_a , 0);
    i = i + 1;
i = 0;

[ I_edges[0][1] , I_edges[0][2] , I_edges[0][3] ]     = attach3( I_edges[0][1] , I_edges[0][2] , I_edges[0][3] )
[ I_edges[2][1] , I_edges[2][2] , I_edges[2][3] ]     = attach3( I_edges[2][1] , I_edges[2][2] , I_edges[2][3] )

[ I_edges[4][1] , I_edges[4][2] , I_edges[4][3] ]     = attach2( I_edges[4][1] , I_edges[4][2] , I_edges[4][3] )


i = 0;
for edge in I_edges:
    edges.append( "  arc  "+I_edges[i][0]+"    ( "+str(I_edges[i][1])+"  "+str(I_edges[i][2])+"  "+str(I_edges[i][3])+" )" );
    i = i + 1;

## Carving
edges.append( "  arc  64 66    ( "+str( c_l+(c_w*0.5) )+"               "+str(MC_l - c_w*(1./2) - c_eD)+"  "+str(MC_h)+" )" );
edges.append( "  arc  65 67    ( "+str( c_l+(c_w*0.5) )+"               "+str(MC_l - c_w*(1./2) - c_eD)+"  "+str(MC_h+c_h)+" )" );
edges.append( "  arc  44 64    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) - c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(MC_h)+" )" );
edges.append( "  arc  45 65    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) - c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(MC_h+c_h)+" )" );
edges.append( "  arc  66 56    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) + c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(MC_h)+" )" );
edges.append( "  arc  67 57    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) + c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(MC_h+c_h)+" )" );

# Carving Sphere
[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , 0);
edges.append( "  arc  57 72    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , -45);
edges.append( "  arc  67 74    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , -45/2.);
edges.append( "  arc  72 74    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , -90);
edges.append( "  arc  74 75    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , -135);
edges.append( "  arc  75 65    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , -180);
edges.append( "  arc  73 45    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45 , -90-45/2.);
edges.append( "  arc  73 75    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( c_l     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   c_l     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 90 , 0);
edges.append( "  arc  72 73    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( 0     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   0     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , 0);
edges.append( "  arc  59 79    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( 0     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   0     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 90 , 0);
edges.append( "  arc  79 82    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

[x,y,z]=rot3d( 0     ,   MC_l - c_w*(1./2) - c_eD    ,   MC_h+c_h   ,   0     ,   MC_l-c_eD   ,   MC_h+c_h  , 0 , 45/2. , -180);
edges.append( "  arc  82 49    ( "+str(x)+"   "+str(y)+"  "+str(z)+" )" );

# NEW main chamber
edges.append( "  arc  87 93    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) - c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(0)+" )" );
edges.append( "  arc  93 94    ( "+str( c_l+(c_w*0.5) )+"               "+str(MC_l - c_w*(1./2) - c_eD)+"  "+str(0)+" )" );
edges.append( "  arc  94 90    ( "+str( c_l+(c_w*0.5)*sinD(45./2) )+"   "+str(MC_l - c_w*(1./2) + c_w*(1./2)*cosD(45./2) - c_eD)+"  "+str(0)+" )" );

edges.append( "  arc  99 97    ( "+str( MC_w-MC_R+MC_R*sinD(45./2) )+"   "+str( MC_l-MC_R+cosD(45/2.)*MC_R )+"  "+str(0)+" )" );
edges.append( "  arc  98 100    ( "+str( MC_w-MC_R+MC_R*sinD(45./2) )+"   "+str( MC_l-MC_R+cosD(45/2.)*MC_R )+"  "+str(MC_h)+" )" );

edges.append( "  spline  95 97 (( "+str( MC_w )+"   "+str( MC_l-MC_R )+"  "+str(0)+" ) ( "+str( MC_w-MC_R+MC_R*cosD(45./2) )+"   "+str( MC_l-MC_R+sinD(45/2.)*MC_R )+"  "+str(0)+" ) ) " ) ;
edges.append( "  spline  96 98 (( "+str( MC_w )+"   "+str( MC_l-MC_R )+"  "+str(MC_h)+" )( "+str( MC_w-MC_R+MC_R*cosD(45./2) )+"   "+str( MC_l-MC_R+sinD(45/2.)*MC_R )+"  "+str(MC_h)+" ) ) " ) ;

# Add points to buffer
# OpenFOAM now rejects obsolete arcs that are not edges of any block.
# These four unused arcs were silently ignored by the original mesher.
edges = [edge for edge in edges if tuple(map(int, edge.split()[1:3]))
         not in {(6,12),(7,13),(12,18),(13,19)}]
for edge in edges:
    buffer = buffer+ edge+"\n";

# Close buffer
buffer = buffer+");\n";

# Write to file
FILE.write(buffer);


###########################################
#       Script - Patches                  #
###########################################
buffer = "";
buffer = """
// Patches, page 38 in notebook

patches
(
    symmetryPlane Symmetri
    (
        // Main Box
        ( 27 58 86 26 )
        ( 58 54 85 86 )
        ( 54 50 84 85 )
        ( 50 48 101 84 )
        ( 48 21 20 101 )

        // Inlet
        ( 35 33 32 34 )
        ( 33 31 30 32 )
        ( 31 29 28 30 )

        // Carving
        ( 51 49 48 50 )
        ( 55 51 50 54 )
        ( 59 55 54 58 )
        ( 55 59 79 80 )
        ( 55 80 81 51 )
        ( 51 81 82 49 )
        ( 80 79 82 81 )

    )

    patch membrane
    (
        ( 101 87 88 84 )
        ( 84 88 89 85 )
        ( 85 89 90 86 )
        ( 87 93 91 88 )
        ( 91 93 94 92 )
        ( 89 92 94 90 )
        ( 88 91 92 89 )
        ( 86 90 99 26 )
        ( 90 94 97 99 )
        ( 93 95 97 94 )
        ( 0 95 93 2 )
        ( 4 2 93 87 )
        ( 20 4 87 101 )
    )


    patch inlet
    (
        ( 31 37 41 29 )
        ( 37 39 43 41 )
        ( 33 35 43 39 )
        ( 31 33 39 37 )
    )

    patch inlet_connector
    (
        // Inlet
        ( 28 40 36 30 )
        ( 40 42 38 36 )
        ( 38 42 34 32 )
        ( 32 30 36 38 )
    )

    wall carvingTop
    (
        ( 59 57 72 79 )
        ( 79 72 73 82 )
        ( 57 59 58 56 )
    )

    wall fixedWalls
    (
        // Top Main chamber
        ( 3 64 96 1 )
        ( 5 44 64 3 )
        ( 21 48 44 5 )
        ( 64 66 98 96 )
        ( 66 56 100 98 )
        ( 56 58 27 100 )

        // Sides
        ( 21 5 4 20 )
        ( 5 3 2 4 )
        ( 3 1 0 2 )
        ( 1 96 95 0 )
        ( 96 98 97 95 )
        ( 98 100 99 97 )
        ( 100 27 26 99 )


        // Carving Walls
        ( 82 73 45 49 )
        ( 49 45 44 48 )
        ( 67 57 72 74 )
        ( 72 73 75 74 )
        ( 73 45 65 75 )
        ( 74 67 65 75 )
        ( 45 65 64 44 )
        ( 65 67 66 64 )
        ( 67 57 56 66 )


        // Inlet
        ( 28 29 41 40 )
        ( 40 41 43 42 )
        ( 42 43 35 34 )

    )
);\n
""";


# Write to file
FILE.write(buffer);

###########################################
#       Script - Merging                  #
###########################################
buffer = "";
buffer = """
// Merging
mergePatchPairs
(
   ( carvingTop inlet_connector )
   //( top carving_connector  )
);
""";

# Write to file
FILE.write(buffer);


###########################################
#       Script - Closing                    #
###########################################
FILE.close();



###########################################
#       Script - OpenFOAM Commands        #
###########################################

# Run blockMesh to create Mesh
