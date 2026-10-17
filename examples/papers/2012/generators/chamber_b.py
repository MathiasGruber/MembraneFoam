# Chamber B dictionary generator, derived from createFO21Dict.py (Mathias Gruber, 2012).
# Run in a new case directory. See ../README.md for usage.
# Import math package
from math import sin, cos, radians, asin, atan, sqrt
from pathlib import Path
Path("system").mkdir(parents=True,exist_ok=True)
import re


###########################################
# Description                             #
# This script creates the blockMeshDict   #
# for creating the FO21 chamber using     #
# blockMesh. Various parameters can be    #
# adjusted using the script. Must be run  #
# from case root directory                #
# All units are in cm                     #
###########################################

###########################################
#       USER INPUT - START                #
###########################################

inlet_diameter = 0.15;
inlet_length = 0.8;
inlet2_length = 0.6;

center_inlet_upper_height = 0.32;
off_center_inlet_upper_height = 0.34;

inlet_vertical_angle = 35.;
inlet_horizontal_angle = 17;

chamber_diameter = 3.0;
chamber_height   = 0.35;


###########################################
#       Advanced input - START            #
###########################################

# Open File
filename = "system/blockMeshDict";
FILE = open(filename,"w");

# Merging Adjustments
inlet_chamber_seperation = 0.001;

# Main chamber adjustments
Mcenter_square_radius = 0.5;

# Center Inlet adjustments
center_square_radius = 0.025;

# Initial center/off-center inlet seperation
iniSep = 0.001;


###########################################
#       Script - HEADER                   #
###########################################
header = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
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

// Dimensions in centi-meters
scale                   0.01;
""";

FILE.write(header);

# Main dimensions in strings
cD = str(chamber_diameter);
cH = str(chamber_height);
cR = str((chamber_diameter / 2));
csR = str(Mcenter_square_radius);
cIU = str(center_inlet_upper_height);

# Main dimensions in floats
fcR = chamber_diameter / 2;
fcIU = center_inlet_upper_height;
foIU = off_center_inlet_upper_height;
fiR = inlet_diameter / 2;

###########################################
#       Script - Functions                #
###########################################

def sinD( x ):
    return sin( radians(x) );

def cosD( x ):
    return cos( radians(x) );

def rotD( oriY , oriZ , poiY , poiZ , angle ):

    yNew = oriY + ( poiY - oriY )* cosD(angle) - ( poiZ - oriZ ) * sinD( angle );
    zNew = oriZ + ( poiY - oriY )* sinD(angle) + ( poiZ - oriZ ) * cosD( angle );
    return [ yNew , zNew ];


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



# Attach central inlet to membrane, return [ y , z ], p.32 in notebook
def attach( posX , posZ , posY ):
    # Determine V
    main_cylinder_Y = cos( asin( posX / fcR ) )*fcR;
    V = main_cylinder_Y - posY;
    # New Positions
    yNew = posY + V + inlet_chamber_seperation;
    zNew = posZ + sinD(inlet_vertical_angle)*V + sinD(inlet_vertical_angle)*inlet_chamber_seperation;
    # New xz positions
    return [ yNew , zNew ];

# Attach Off-central inlet to membrane, return [ x, y , z ], p.33_left in notebook
def attach2( posX ,  posY , posZ):
    # Determine move vector
    vector = rot3d( 0 , 0 , 0 , 0 , inlet_chamber_seperation , 0 , 0 ,inlet_vertical_angle , inlet_horizontal_angle );
    #print vector;

    # Determine distance
    distance = inlet_chamber_seperation*2;
    eps = inlet_chamber_seperation; i = 0;
    while ( eps < distance and i < 1000):
        posX = posX - vector[0]
        posY = posY - vector[1]
        posZ = posZ - vector[2]
        # Determine distance
        alpha = atan( posX / posY );
        circleX_Y = [ sin(alpha)*(fcR+inlet_chamber_seperation) , cos(alpha)*(fcR+inlet_chamber_seperation) ]
        distance = sqrt((posX - circleX_Y[0])**2 + (posY - circleX_Y[1])**2)
        #print distance;

    # New xz positions
    return [ posX, posY , posZ ];

###########################################
#       Script - Vertices                 #
###########################################

print("\nCalculating vertice positions for: Main Cylinder ");

# Numbers used underway

# Main Cylinder
xy_CWall = str(sinD(45)*fcR);

print("\nCalculating vertice positions for: Center Inlet ");

# Center Inlet
x_15_19_22_24 = sinD(45)*fiR;
x_16_20_23_25 = center_square_radius;

z_14          = fcIU;
z_18          = fcIU;
z_15          = fcIU-(fiR-cosD(45)*fiR) ;
z_19          = fcIU-(fiR-cosD(45)*fiR) ;
z_22          = fcIU - fiR - cosD(45)*fiR ;
z_24          = fcIU - fiR - cosD(45)*fiR ;
z_17_16       = fcIU - (fiR - center_square_radius) ;
z_21_20       = fcIU - (fiR - center_square_radius) ;
z_27_23       = fcIU - (fiR + center_square_radius) ;
z_29_25       = fcIU - (fiR + center_square_radius) ;
z_26          = fcIU - inlet_diameter
z_28          = fcIU - inlet_diameter

y_14       = chamber_diameter / 2 + inlet_chamber_seperation
y_15       = chamber_diameter / 2 + inlet_chamber_seperation
y_17_16    = chamber_diameter / 2 + inlet_chamber_seperation

y_27_23    = chamber_diameter / 2 + inlet_chamber_seperation
y_22       = chamber_diameter / 2 + inlet_chamber_seperation
y_26       = chamber_diameter / 2 + inlet_chamber_seperation

y_18       = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
y_19       = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
y_21_20    = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation

y_29_25    = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
y_28       = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
y_24       = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation

# Points of rotation & rotations
cIrY = chamber_diameter / 2 + inlet_chamber_seperation;
cIrZ = center_inlet_upper_height - fiR;

[ y_14 , z_14]          = rotD( cIrY , cIrZ , y_14 , z_14 , inlet_vertical_angle );
[ y_15 , z_15]          = rotD( cIrY , cIrZ , y_15 , z_15 , inlet_vertical_angle );
[ y_17_16 , z_17_16]    = rotD( cIrY , cIrZ , y_17_16 , z_17_16 , inlet_vertical_angle );

[ y_27_23 , z_27_23]    = rotD( cIrY , cIrZ , y_27_23 , z_27_23 , inlet_vertical_angle );
[ y_22 , z_22]          = rotD( cIrY , cIrZ , y_22 , z_22 , inlet_vertical_angle );
[ y_26 , z_26]          = rotD( cIrY , cIrZ , y_26 , z_26 , inlet_vertical_angle );

[ y_18 , z_18]          = rotD( cIrY , cIrZ , y_18 , z_18 , inlet_vertical_angle );
[ y_19 , z_19]          = rotD( cIrY , cIrZ , y_19 , z_19 , inlet_vertical_angle );
[ y_21_20 , z_21_20]    = rotD( cIrY , cIrZ , y_21_20, z_21_20 , inlet_vertical_angle );

[ y_29_25 , z_29_25]    = rotD( cIrY , cIrZ , y_29_25, z_29_25 , inlet_vertical_angle );
[ y_28 , z_28]          = rotD( cIrY , cIrZ , y_28, z_28 , inlet_vertical_angle );
[ y_24 , z_24]          = rotD( cIrY , cIrZ , y_24, z_24 , inlet_vertical_angle );

# Attaching to flow chamber
[ y_14 , z_14]          = attach( 0             , z_14      , y_14 );
[ y_15 , z_15]          = attach( x_15_19_22_24 , z_15      , y_15 );
[ y_17_16 , z_17_16]    = attach( fiR/2         , z_17_16   , y_17_16 );

[ y_27_23 , z_27_23]    = attach( fiR/2         , z_27_23   , y_27_23 );
[ y_22 , z_22]          = attach( x_15_19_22_24 , z_22      , y_22 );
[ y_26 , z_26]          = attach( 0             , z_26      , y_26  );

print("\nCalculating vertice positions for: Off-center Inlet ");
# Off-Center Inlet
x_30 = fiR + ( fiR - center_square_radius ) + iniSep;
x_31 = fiR + ( fiR + center_square_radius ) + iniSep;
x_41 = fiR + ( fiR - cosD(45)*fiR ) + iniSep;
x_38 = fiR + ( fiR + cosD(45)*fiR ) + iniSep;
x_33 = fiR + ( fiR - center_square_radius ) + iniSep;
x_32 = fiR + ( fiR + center_square_radius ) + iniSep;
x_40 = fiR + ( fiR - cosD(45)*fiR ) + iniSep;
x_39 = fiR + ( fiR + cosD(45)*fiR ) + iniSep;

x_34 = fiR + ( fiR - center_square_radius ) + iniSep;
x_35 = fiR + ( fiR + center_square_radius ) + iniSep;
x_45 = fiR + ( fiR - cosD(45)*fiR ) + iniSep;
x_42 = fiR + ( fiR + cosD(45)*fiR ) + iniSep;
x_37 = fiR + ( fiR - center_square_radius ) + iniSep;
x_36 = fiR + ( fiR + center_square_radius ) + iniSep;
x_44 = fiR + ( fiR - cosD(45)*fiR ) + iniSep;
x_43 = fiR + ( fiR + cosD(45)*fiR ) + iniSep;

z_41 = foIU - ( fiR - sinD(45)*fiR );
z_30 = foIU - ( fiR - center_square_radius );
z_33 = foIU - ( fiR + center_square_radius );
z_40 = foIU - ( fiR + sinD(45)*fiR )
z_38 = foIU - ( fiR - sinD(45)*fiR );
z_31 = foIU - ( fiR - center_square_radius );
z_32 = foIU - ( fiR + center_square_radius );
z_39 = foIU - ( fiR + sinD(45)*fiR )

z_45 = foIU - ( fiR - sinD(45)*fiR );
z_34 = foIU - ( fiR - center_square_radius );
z_37 = foIU - ( fiR + center_square_radius );
z_44 = foIU - ( fiR + sinD(45)*fiR );
z_42 = foIU - ( fiR - sinD(45)*fiR );
z_35 = foIU - ( fiR - center_square_radius );
z_36 = foIU - ( fiR + center_square_radius );
z_43 = foIU - ( fiR + sinD(45)*fiR );

y_30 = chamber_diameter / 2 + inlet_chamber_seperation;
y_31 = chamber_diameter / 2 + inlet_chamber_seperation;
y_32 = chamber_diameter / 2 + inlet_chamber_seperation;
y_33 = chamber_diameter / 2 + inlet_chamber_seperation;
y_38 = chamber_diameter / 2 + inlet_chamber_seperation;
y_39 = chamber_diameter / 2 + inlet_chamber_seperation;
y_40 = chamber_diameter / 2 + inlet_chamber_seperation;
y_41 = chamber_diameter / 2 + inlet_chamber_seperation;


y_34 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_35 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_36 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_37 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_42 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_43 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_44 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
y_45 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;

offCenterInletVertices = [];
offCenterInletVertices.append( [ x_30 , y_30 , z_30 ] );
offCenterInletVertices.append( [ x_31 , y_31 , z_31 ] );
offCenterInletVertices.append( [ x_32 , y_32 , z_32 ] );
offCenterInletVertices.append( [ x_33 , y_33 , z_33 ] );
offCenterInletVertices.append( [ x_34 , y_34 , z_34 ] );
offCenterInletVertices.append( [ x_35 , y_35 , z_35 ] );
offCenterInletVertices.append( [ x_36 , y_36 , z_36 ] );
offCenterInletVertices.append( [ x_37 , y_37 , z_37 ] );
offCenterInletVertices.append( [ x_38 , y_38 , z_38 ] );
offCenterInletVertices.append( [ x_39 , y_39 , z_39 ] );
offCenterInletVertices.append( [ x_40 , y_40 , z_40 ] );
offCenterInletVertices.append( [ x_41 , y_41 , z_41 ] );
offCenterInletVertices.append( [ x_42 , y_42 , z_42 ] );
offCenterInletVertices.append( [ x_43 , y_43 , z_43 ] );
offCenterInletVertices.append( [ x_44 , y_44 , z_44 ] );
offCenterInletVertices.append( [ x_45 , y_45 , z_45 ] );


# points of rotation & rotation: yz plane
ocIrY = chamber_diameter / 2 + inlet_chamber_seperation;
ocIrZ = off_center_inlet_upper_height - fiR;

i = 0;
for point in offCenterInletVertices:
    offCenterInletVertices[i] = rot3d( 0 , ocIrY , ocIrZ , offCenterInletVertices[i][0],  offCenterInletVertices[i][1] ,  offCenterInletVertices[i][2] , 0 , inlet_vertical_angle , 0);
    i = i + 1;


# points of rotation & rotation: yx plane
ocIr2X = fiR + iniSep;
ocIr2Y = fcR + cosD( inlet_vertical_angle )*inlet2_length;

i = 0;
for point in offCenterInletVertices:
    offCenterInletVertices[i] = rot3d( ocIr2X , ocIr2Y , 0 , offCenterInletVertices[i][0],  offCenterInletVertices[i][1] ,  offCenterInletVertices[i][2] , 0 , 0 , inlet_horizontal_angle);
    i = i + 1;

# Attach to membrane
offCenterInletVertices[0] = attach2( offCenterInletVertices[0][0],  offCenterInletVertices[0][1] ,  offCenterInletVertices[0][2] );
offCenterInletVertices[1] = attach2( offCenterInletVertices[1][0],  offCenterInletVertices[1][1] ,  offCenterInletVertices[1][2] );
offCenterInletVertices[2] = attach2( offCenterInletVertices[2][0],  offCenterInletVertices[2][1] ,  offCenterInletVertices[2][2] );
offCenterInletVertices[3] = attach2( offCenterInletVertices[3][0],  offCenterInletVertices[3][1] ,  offCenterInletVertices[3][2] );
offCenterInletVertices[8] = attach2( offCenterInletVertices[8][0],  offCenterInletVertices[8][1] ,  offCenterInletVertices[8][2] );
offCenterInletVertices[9] = attach2( offCenterInletVertices[9][0],  offCenterInletVertices[9][1] ,  offCenterInletVertices[9][2] );
offCenterInletVertices[10] = attach2( offCenterInletVertices[10][0],  offCenterInletVertices[10][1] ,  offCenterInletVertices[10][2] );
offCenterInletVertices[11] = attach2( offCenterInletVertices[11][0],  offCenterInletVertices[11][1] ,  offCenterInletVertices[11][2] );


buffer = """
// Points in {xyz}
vertices
(
""";

# MAIN CHAMBER POINTS
vertices = [];
vertices.append( "    ( 0 0 0)                              // Number 0" );
vertices.append( "    ( "+csR+" 0 0 )                       // Number 1" );
vertices.append( "    ( "+csR+" "+csR+" 0 )                 // Number 2" );
vertices.append( "    ( 0 "+csR+" 0 )                       // Number 3" );
vertices.append( "    ( 0 0 "+cH+")                         // Number 4" );
vertices.append( "    ( "+csR+" 0 "+cH+" )                  // Number 5" );
vertices.append( "    ( "+csR+" "+csR+" "+cH+" )            // Number 6" );
vertices.append( "    ( 0 "+csR+" "+cH+" )                  // Number 7" );
vertices.append( "    ( "+cR+" 0 0 )                        // Number 8" );
vertices.append( "    ( "+xy_CWall+" "+xy_CWall+" 0 )       // Number 9" );
vertices.append( "    ( "+cR+" 0 "+cH+" )                   // Number 10" );
vertices.append( "    ( "+xy_CWall+" "+xy_CWall+" "+cH+" )  // Number 11" );
vertices.append( "    ( 0 "+cR+" 0 )                        // Number 12" );
vertices.append( "    ( 0 "+cR+" "+cH+" )                   // Number 13" );

# Center inlet
vertices.append( "    ( 0                           "+str(y_14)+"            "+str(z_14)+" )              // Number 14" );
vertices.append( "    ( "+str(x_15_19_22_24)+"      "+str(y_15)+"            "+str(z_15)+" )              // Number 15" );
vertices.append( "    ( "+str(x_16_20_23_25)+"      "+str(y_17_16)+"         "+str(z_17_16)+" )           // Number 16" );
vertices.append( "    ( 0                           "+str(y_17_16)+"         "+str(z_17_16)+" )           // Number 17" );
vertices.append( "    ( 0                           "+str(y_18)+"            "+str(z_18)+" )              // Number 18" );
vertices.append( "    ( "+str(x_15_19_22_24)+"      "+str(y_19)+"            "+str(z_19)+" )              // Number 19" );
vertices.append( "    ( "+str(x_16_20_23_25)+"      "+str(y_21_20)+"         "+str(z_21_20)+" )           // Number 20" );
vertices.append( "    ( 0                           "+str(y_21_20)+"         "+str(z_21_20)+" )           // Number 21" );

vertices.append( "    ( "+str(x_15_19_22_24)+"      "+str(y_22)+"            "+str(z_22)+")                // Number 22" );
vertices.append( "    ( "+str(x_16_20_23_25)+"      "+str(y_27_23)+"         "+str(z_27_23)+"    )         // Number 23" );
vertices.append( "    ( "+str(x_15_19_22_24)+"      "+str(y_24)+"            "+str(z_24)+"    )            // Number 24" );
vertices.append( "    ( "+str(x_16_20_23_25)+"      "+str(y_29_25)+"         "+str(z_29_25)+"    )         // Number 25" );
vertices.append( "    ( 0                           "+str(y_26)+"            "+str(z_26)+"    )            // Number 26" );
vertices.append( "    ( 0                           "+str(y_27_23)+"         "+str(z_27_23)+"    )         // Number 27" );
vertices.append( "    ( 0                           "+str(y_28)+"            "+str(z_28)+"    )            // Number 28" );
vertices.append( "    ( 0                           "+str(y_29_25)+"         "+str(z_29_25)+"    )         // Number 29" );

# Off-Center inlet
i = 30;
for vertex in offCenterInletVertices:
    vertices.append( "    ( "+str(vertex[0])+"    "+str(vertex[1])+"    "+str(vertex[2])+"    )   // Number "+str(i) );
    i = i + 1;



# Add points to buffer
for vertice in vertices:
    buffer = buffer+ vertice+"\n";

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

# Main chamber edges
edges.append( "    arc 8 9      ( "+str(cosD(45/2)*fcR)+" "+str(sinD(45/2)*fcR)+" 0 )" );
edges.append( "    arc 10 11    ( "+str(cosD(45/2)*fcR)+" "+str(sinD(45/2)*fcR)+" "+cH+" )" );
edges.append( "    arc 9 12     ( "+str(sinD(45/2)*fcR)+" "+str(cosD(45/2)*fcR)+" 0 )" );
edges.append( "    arc 11 13    ( "+str(sinD(45/2)*fcR)+" "+str(cosD(45/2)*fcR)+" "+cH+" )" );

# Center Inlet edges
XedgeP_middle = sinD(22.5)*fiR;
XedgeP_15_19_22_24             = fiR;

YedgeP_14_15 = chamber_diameter / 2 + inlet_chamber_seperation
YedgeP_18_19 = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
YedgeP_15_22 = chamber_diameter / 2 + inlet_chamber_seperation
YedgeP_19_24 = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation
YedgeP_22_26 = chamber_diameter / 2 + inlet_chamber_seperation
YedgeP_24_28 = chamber_diameter / 2 + inlet_length + inlet_chamber_seperation

ZedgeP_14_15 = fcIU - ( fiR - cosD(22.5)*fiR)
ZedgeP_18_19 = fcIU - ( fiR - cosD(22.5)*fiR)
ZedgeP_15_22 = fcIU - fiR
ZedgeP_19_24 = fcIU - fiR
ZedgeP_22_26 = fcIU - fiR - cosD(22.5)*fiR
ZedgeP_24_28 = fcIU - fiR - cosD(22.5)*fiR

# Apply Rotations
[ YedgeP_14_15 , ZedgeP_14_15]          = rotD( cIrY , cIrZ , YedgeP_14_15 , ZedgeP_14_15 , inlet_vertical_angle );
[ YedgeP_18_19 , ZedgeP_18_19]          = rotD( cIrY , cIrZ , YedgeP_18_19 , ZedgeP_18_19 , inlet_vertical_angle );
[ YedgeP_15_22 , ZedgeP_15_22]          = rotD( cIrY , cIrZ , YedgeP_15_22 , ZedgeP_15_22 , inlet_vertical_angle );
[ YedgeP_19_24 , ZedgeP_19_24]          = rotD( cIrY , cIrZ , YedgeP_19_24 , ZedgeP_19_24 , inlet_vertical_angle );
[ YedgeP_22_26 , ZedgeP_22_26]          = rotD( cIrY , cIrZ , YedgeP_22_26 , ZedgeP_22_26 , inlet_vertical_angle );
[ YedgeP_24_28 , ZedgeP_24_28]          = rotD( cIrY , cIrZ , YedgeP_24_28 , ZedgeP_24_28 , inlet_vertical_angle );

# Attach to chamber
[ YedgeP_14_15 , ZedgeP_14_15]          = attach( XedgeP_middle         , ZedgeP_14_15      , YedgeP_14_15 );
[ YedgeP_15_22 , ZedgeP_15_22]          = attach( XedgeP_15_19_22_24    , ZedgeP_15_22      , YedgeP_15_22 );
[ YedgeP_22_26 , ZedgeP_22_26]          = attach( XedgeP_middle         , ZedgeP_22_26      , YedgeP_22_26 );

edges.append( "    arc 14 15   (   "+str(XedgeP_middle)+"            "+str(YedgeP_14_15)+"                 "+str(ZedgeP_14_15)+"    ) " );
edges.append( "    arc 18 19   (   "+str(XedgeP_middle)+"            "+str(YedgeP_18_19)+"                 "+str(ZedgeP_18_19)+"    ) " );
edges.append( "    arc 15 22   (   "+str(XedgeP_15_19_22_24)+"       "+str(YedgeP_15_22)+"                 "+str(ZedgeP_15_22)+"    ) " );
edges.append( "    arc 19 24   (   "+str(XedgeP_15_19_22_24)+"       "+str(YedgeP_19_24)+"                 "+str(ZedgeP_19_24)+"    ) " );
edges.append( "    arc 22 26   (   "+str(XedgeP_middle)+"            "+str(YedgeP_22_26)+"                 "+str(ZedgeP_22_26)+"    ) " );
edges.append( "    arc 24 28   (   "+str(XedgeP_middle)+"            "+str(YedgeP_24_28)+"                 "+str(ZedgeP_24_28)+"    ) " );

# Off-center inlet edges, initial
Xedge38_39 = 3*fiR + iniSep;
Xedge42_43 = 3*fiR + iniSep;
Xedge39_40 = 2*fiR + iniSep;
Xedge43_44 = 2*fiR + iniSep;
Xedge40_41 = fiR + iniSep;
Xedge44_45 = fiR + iniSep;
Xedge41_38 = 2*fiR + iniSep;
Xedge45_42 = 2*fiR + iniSep;

Yedge38_39 = chamber_diameter / 2 + inlet_chamber_seperation;
Yedge42_43 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
Yedge39_40 = chamber_diameter / 2 + inlet_chamber_seperation;
Yedge43_44 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
Yedge40_41 = chamber_diameter / 2 + inlet_chamber_seperation;
Yedge44_45 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;
Yedge41_38 = chamber_diameter / 2 + inlet_chamber_seperation;
Yedge45_42 = chamber_diameter / 2 + inlet2_length + inlet_chamber_seperation;

Zedge38_39 = foIU - fiR;
Zedge42_43 = foIU - fiR;
Zedge39_40 = foIU - 2*fiR;
Zedge43_44 = foIU - 2*fiR;
Zedge40_41 = foIU - fiR;
Zedge44_45 = foIU - fiR;
Zedge41_38 = foIU;
Zedge45_42 = foIU;

# Apply Rotations: yz
[ Xedge38_39 , Yedge38_39 , Zedge38_39 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge38_39,  Yedge38_39 ,  Zedge38_39 , 0 , inlet_vertical_angle , 0);
[ Xedge42_43 , Yedge42_43 , Zedge42_43 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge42_43,  Yedge42_43 ,  Zedge42_43 , 0 , inlet_vertical_angle , 0);
[ Xedge39_40 , Yedge39_40 , Zedge39_40 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge39_40,  Yedge39_40 ,  Zedge39_40 , 0 , inlet_vertical_angle , 0);
[ Xedge43_44 , Yedge43_44 , Zedge43_44 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge43_44,  Yedge43_44 ,  Zedge43_44 , 0 , inlet_vertical_angle , 0);
[ Xedge40_41 , Yedge40_41 , Zedge40_41 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge40_41,  Yedge40_41 ,  Zedge40_41 , 0 , inlet_vertical_angle , 0);
[ Xedge44_45 , Yedge44_45 , Zedge44_45 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge44_45,  Yedge44_45 ,  Zedge44_45 , 0 , inlet_vertical_angle , 0);
[ Xedge41_38 , Yedge41_38 , Zedge41_38 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge41_38,  Yedge41_38 ,  Zedge41_38 , 0 , inlet_vertical_angle , 0);
[ Xedge45_42 , Yedge45_42 , Zedge45_42 ] = rot3d( 0 , ocIrY , ocIrZ , Xedge45_42,  Yedge45_42 ,  Zedge45_42 , 0 , inlet_vertical_angle , 0);

[ Xedge38_39 , Yedge38_39 , Zedge38_39 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge38_39,  Yedge38_39 ,  Zedge38_39 , 0 , 0 , inlet_horizontal_angle);
[ Xedge42_43 , Yedge42_43 , Zedge42_43 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge42_43,  Yedge42_43 ,  Zedge42_43 , 0 , 0 , inlet_horizontal_angle);
[ Xedge39_40 , Yedge39_40 , Zedge39_40 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge39_40,  Yedge39_40 ,  Zedge39_40 , 0 , 0 , inlet_horizontal_angle);
[ Xedge43_44 , Yedge43_44 , Zedge43_44 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge43_44,  Yedge43_44 ,  Zedge43_44 , 0 , 0 , inlet_horizontal_angle);
[ Xedge40_41 , Yedge40_41 , Zedge40_41 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge40_41,  Yedge40_41 ,  Zedge40_41 , 0 , 0 , inlet_horizontal_angle);
[ Xedge44_45 , Yedge44_45 , Zedge44_45 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge44_45,  Yedge44_45 ,  Zedge44_45 , 0 , 0 , inlet_horizontal_angle);
[ Xedge41_38 , Yedge41_38 , Zedge41_38 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge41_38,  Yedge41_38 ,  Zedge41_38 , 0 , 0 , inlet_horizontal_angle);
[ Xedge45_42 , Yedge45_42 , Zedge45_42 ] = rot3d( ocIr2X , ocIr2Y , 0 , Xedge45_42,  Yedge45_42 ,  Zedge45_42 , 0 , 0 , inlet_horizontal_angle);

[ Xedge38_39 , Yedge38_39 , Zedge38_39] = attach2( Xedge38_39         , Yedge38_39      , Zedge38_39 );
[ Xedge40_41 , Yedge40_41 , Zedge40_41] = attach2( Xedge40_41         , Yedge40_41      , Zedge40_41 );
[ Xedge39_40 , Yedge39_40 , Zedge39_40] = attach2( Xedge39_40         , Yedge39_40      , Zedge39_40 );
[ Xedge41_38 , Yedge41_38 , Zedge41_38] = attach2( Xedge41_38         , Yedge41_38      , Zedge41_38 );

edges.append( "    arc 38 39    ( "+str(Xedge38_39)+" "+str(Yedge38_39)+" "+str(Zedge38_39)+" )" );
edges.append( "    arc 42 43    ( "+str(Xedge42_43)+" "+str(Yedge42_43)+" "+str(Zedge42_43)+" )" );
edges.append( "    arc 39 40    ( "+str(Xedge39_40)+" "+str(Yedge39_40)+" "+str(Zedge39_40)+" )" );
edges.append( "    arc 43 44    ( "+str(Xedge43_44)+" "+str(Yedge43_44)+" "+str(Zedge43_44)+" )" );
edges.append( "    arc 40 41    ( "+str(Xedge40_41)+" "+str(Yedge40_41)+" "+str(Zedge40_41)+" )" );
edges.append( "    arc 44 45    ( "+str(Xedge44_45)+" "+str(Yedge44_45)+" "+str(Zedge44_45)+" )" );
edges.append( "    arc 41 38    ( "+str(Xedge41_38)+" "+str(Yedge41_38)+" "+str(Zedge41_38)+" )" );
edges.append( "    arc 45 42    ( "+str(Xedge45_42)+" "+str(Yedge45_42)+" "+str(Zedge45_42)+" )" );



# Add points to buffer
for edge in edges:
    buffer = buffer+ edge+"\n";

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

# Main Chamber Block
blocks.append( "    hex ( 0 1 2 3 4 5 6 7 ) ( 25 25 70 ) simpleGrading ( 1 1 5 )" );
blocks.append( "    hex ( 1 8 9 2 5 10 11 6 ) ( 25 25 70 ) simpleGrading ( 1 1 5 )" );
blocks.append( "    hex ( 2 9 12 3 6 11 13 7 ) ( 25 25 70 ) simpleGrading ( 1 1 5 )" );

# Center Inlet
blocks.append( "    hex ( 14 15 16 17 18 19 20 21 ) ( 3 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 16 15 22 23 20 19 24 25 ) ( 7 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 27 23 22 26 29 25 24 28 ) ( 3 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 17 16 23 27 21 20 25 29 ) ( 3 7 20 ) simpleGrading (1 1 1)" );

# Off-center Inlet
blocks.append( "    hex ( 30 31 32 33 34 35 36 37 ) ( 7 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 30 41 38 31 34 45 42 35 ) ( 7 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 38 39 32 31 42 43 36 35 ) ( 7 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 32 39 40 33 36 43 44 37 ) ( 7 7 20 ) simpleGrading (1 1 1)" );
blocks.append( "    hex ( 40 41 30 33 44 45 34 37 ) ( 7 7 20 ) simpleGrading (1 1 1)" );

# Add points to buffer
for block in blocks:
    buffer = buffer+ block+"\n";

# Close buffer
buffer = buffer+");\n";

# Write to file
FILE.write(buffer);

###########################################
#       Script - Patches                  #
###########################################
buffer = "";
buffer = """
// Patches
patches
(
""";

blocks = [];

# Main Chamber Block
blocks.append( """
    // Connection with Inlet
    wall cylinderInlets
    (
        ( 9 12 13 11 )
    )

    patch inlet
    (
        // Center Inlet
        ( 18 19 20 21 )
        ( 19 24 25 20 )
        ( 24 28 29 25 )
        ( 29 21 20 25 )

        // Off-center Inlet
        ( 45 42 35 34 )
        ( 42 43 36 35 )
        ( 43 44 37 36 )
        ( 44 45 34 37 )
        ( 34 35 36 37 )
    )

    symmetryPlane Symmetri
    (
        ( 7 4 0 3 )
        ( 7 3 12 13 )
        ( 26 27 29 28 )
        ( 27 17 21 29 )
        ( 17 14 18 21 )
    )

    wall inletConnection
    (
       // Center Cylinder
       ( 17 16 15 14 )
       ( 22 15 16 23 )
       ( 26 22 23 27 )
       ( 27 23 16 17 )

       // Off-center Cylinder
       ( 38 41 30 31 )
       ( 31 32 39 38 )
       ( 32 33 40 39 )
       ( 40 33 30 41 )
       ( 30 33 32 31 )
    )

    wall fixedWalls
        (
            // Cylinder Walls
            ( 8 9 11 10 )

            // Internal Wall
            ( 10 5 1 8 )
            ( 5 4 0 1 )

            // Top
            ( 5 6 7 4 )
            ( 10 11 6 5 )
            (6 11 13 7 )


            // Center Cylinder
            ( 14 15 19 18 )
            ( 15 22 24 19 )
            ( 22 26 28 24 )

            // Off Center Cylinder
            ( 41 38 42 45 )
            ( 38 39 43 42 )
            ( 39 40 44 43 )
            ( 40 41 45 44 )
        )

    patch wallMembrane
        (
            // Membrane
            ( 0 3 2 1 )
            ( 1 2 9 8 )
            ( 2 3 12 9 )
        )
""" );


# Add points to buffer
for block in blocks:
    buffer = buffer+ block+"\n";

# Close buffer
buffer = buffer+");\n";

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
""";

merges = [];

# Main Chamber Block
merges.append( "( inletConnection cylinderInlets )" );


# Add points to buffer
for merge in merges:
    buffer = buffer+ merge+"\n";

# Close buffer
buffer = buffer+");\n";

# Write to file
FILE.write(buffer);

# Close the file being written
FILE.close();


###########################################
#       Script - OpenFOAM Commands        #
###########################################

# Run blockMesh to create Mesh
