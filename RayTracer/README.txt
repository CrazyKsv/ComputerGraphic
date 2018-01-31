Notes:

 * Don't forget to complete the autotester_id.txt
 * If you have done any crunchy work, put your rendered images (in .jpg format, good resolution)
   along with a text file describing your crunchy features in a .zip file called 'crunchy.zip'
   and submit that along with the remaning requested files and/or renders.

What have been done so far:

1. basic texture mapping, Antialiasing, Area light sources (all theses feature support Sphere, Cylinder and Plane)
   supported multi-areaLightSources
2.Normal map. code file : utils.c line 326-335-plane, 401-425-sphere and 492-516-Cylinder
3.Alpha map. code file :  utils.c line 899-924
4.Hierarchical object. Code file : RayTracer.c line 60-83. Support sphere and Cylinder
5.Photon map. included tracePhotonRay,setPixelCol,boxBlur functions and main code from line 638-680
6.Multi-thread change in compile.sh, 
  add -fopenmp and 
  #pragma omp parallel for schedule(dynamic, 32) shared(anti_k, rgbIm, object_list) private(j)
  #pragma omp parallel for private(pc, d, ray, col, new_col, i)
  before main loop.

Feature image:
hierarchical_cyl.ppm is hierarchical object with respect to Cylinder (tree structure)
hierarchical_sphere.ppm is hierarchical object with respect to Sphere (tree structure)
hierarchical_with_refraction.pmm is single structure with refracted sphere

refraction.ppm has a sphere/Ball that is refracted and looks refracted.

normal_alpha_map_photon has normal map on both sphere and plane.
                        has alpha map on the back plane.
                        has photon map on the shadow
normal_alpha_map_photon2 has normal map on two planes.
                        has alpha map on the bottom plane.
                        has photon map on the shadow
photon.ppm is the photon mapping.
twoLightSources.ppm has two area light source with different color light.

finalFancyChrunchy: has texture and normal mapping in the two air shperes and refraction for sphere on the floor
Photon map shows behind the ball in the shadow can be seen in mirrow planes.
The greatest man paco stares at the whole universe. This must be the most fancy moment. Put in the centre between two back planes.
We have a normal mapped photo frame for paco.
2 Mirrow plane with global value 1. Back and top/bottom planes with global value 0.
The ball on the floor has ALPHA value 0.
Light source is close to the camera and is behind the camera.
4 planes consists of a rhombus.




 
