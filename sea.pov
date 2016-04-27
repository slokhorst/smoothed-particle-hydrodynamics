// PoVRay 3.7 Scene File "Waterglass_01c.pov"
// author:  Friedrich A. Lohmueller, March-2011
// Homepage: http://www.f-lohmueller.de
//------------------------------------------------------------------------
#version 3.7;
global_settings{ assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9 }} 
//------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"

#debug concat("input: data/posdata",str(frame_number,-3,0),".csv\n")

//------------------------------------------------------------------------
#declare Camera_0 = camera {/*ultra_wide_angle*/ angle 15      // front view
                            location  <0.0 , 1.0 ,-40.0>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}
#declare Camera_1 = camera {/*ultra_wide_angle*/ angle 90 // diagonal view
                            location  <-2.0 , 1.0 , -2.0>
                            right     x*image_width/image_height
                            look_at   < 0.0 , 1.0 , 0.0>}
#declare Camera_2 = camera {/*ultra_wide_angle*/ angle 90  //right side view
                            location  <3.0 , 1.0 , 0.0>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}
#declare Camera_3 = camera {/*ultra_wide_angle*/ angle 90        // top view
                            location  <0.0 , 3.0 ,-0.001>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}
camera{Camera_1}
//------------------------------------------------------------------------
// sun -------------------------------------------------------------------
light_source{<-1500,2500,-2500> color White*0.9 }
// sky -------------------------------------------------------------------
sky_sphere{ pigment{ gradient <0,1,0>
                     color_map{ [0   color rgb<1,1,1>         ]//White
                                [0.4 color rgb<0.24,0.34,0.56>*0.8]//~Navy
                                [0.6 color rgb<0.24,0.34,0.56>*0.8]//~Navy
                                [1.0 color rgb<1,1,1>         ]//White
                              }
                     scale 2 }
           } // end of sky_sphere 
//------------------------------------------------------------------------

// ground -----------------------------------------------------------------
//---------------------------------<<< settings of squared plane dimensions
#declare RasterScale = 0.90;
#declare RasterHalfLine  = 0.0125;  
#declare RasterHalfLineZ = 0.0125; 
//-------------------------------------------------------------------------
#macro Raster(RScale, HLine) 
       pigment{ gradient x scale RScale
                color_map{[0.000   color rgbt<1,1,1,0>*0.8]
                          [0+HLine color rgbt<1,1,1,0>*0.8]
                          [0+HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,0>*0.8]
                          [1.000   color rgbt<1,1,1,0>*0.8]} }
 #end// of Raster(RScale, HLine)-macro    
//-------------------------------------------------------------------------
    
// squared plane XZ
plane { <0,1,0>, 0    // plane with layered textures
        texture { pigment{color White*1.2}
                }
        texture { Raster(RasterScale,RasterHalfLine ) rotate<0,0,0> }
        texture { Raster(RasterScale,RasterHalfLineZ) rotate<0,90,0>}
        rotate<0,0,0>
      }
//------------------------------------------------ end of squared plane XZ

//--------------------------------------------------------------------------
//---------------------------- objects in scene ----------------------------
//--------------------------------------------------------------------------

union{ 
  #declare Val1 = 0.0;
  #declare Val2 = 0.0;
  #declare Val3 = 0.0;
  #fopen MyFile concat("data/posdata",str(frame_number,-3,0),".csv") read
  #while (defined(MyFile))
    #read (MyFile,Val1,Val3,Val2)

    sphere{<0.1*Val1,0.1*(Val2+10),0.1*Val3>,0.08  
        material{
         texture{
          pigment{ rgbf<.93,.95,.98,0.9>*0.95}
          normal { ripples 0.135 scale 0.0125 turbulence 0.01 translate<-0.05,0,0> rotate<0,-0.20,0>} 
          finish { ambient 0.0
                   diffuse 0.15
                   reflection 0.2
                   specular 0.6
                   roughness 0.005
                  // phong 1 
                  // phong_size 400
                   reflection { 0.17, 1.0 fresnel on }
                   conserve_energy
                 }
           } // end of texture
          interior{ ior 1.33 
                    fade_power 1001
                    fade_distance 0.5
                    fade_color <0.2,0.2,0.9> 
                } // end of interior
        } // end of material
    } //
      
 #end // ----------- end of #for loop

} //------------------------------------------------------

