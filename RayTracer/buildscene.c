// Sets up all objects in the scene. This involves creating each object,
// defining the transformations needed to shape and position it as
// desired, specifying the reflectance properties (albedos and colours)
// and setting up textures where needed.
// Light sources must be defined, positioned, and their colour defined.
// All objects must be inserted in the object_list. All light sources
// must be inserted in the light_list.
//
// To create hierarchical objects:
//    You must keep track of transformations carried out by parent objects
//    as you move through the hierarchy. Declare and manipulate your own
//    transformation matrices (use the provided functions in utils.c to
//    compound transformations on these matrices). When declaring a new
//    object within the hierarchy
//    - Initialize the object
//    - Apply any object-level transforms to shape/rotate/resize/move
//      the object using regular object transformation functions
//    - Apply the transformations passed on from the parent object
//      by pre-multiplying the matrix containing the parent's transforms
//      with the object's own transformation matrix.
//    - Compute and store the object's inverse transform as usual.
//
// NOTE: After setting up the transformations for each object, don't
//       forget to set up the inverse transform matrix!

struct object3D *o;
struct pointLS *l;
struct point3D p;

// Simple scene for Assignment 3:
// Insert a couple of objects. A plane and two spheres
// with some transformations.

// Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

//  o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);		// Initialize a sphere
//  Scale(o,1.5,.75,.75);					// Apply a few transforms (Translate * Rotate * Scale)
//  RotateZ(o,PI/4);					
//  Translate(o,2.0,2.5,1.5);
// //  loadTexture(o,"./Texture/ocean.ppm",1,&texture_list);  
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  // If needed, this is how you load a texture map
//  // loadTexture(o,"./Texture/mosaic2.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
// 								// texture gets added to the texture list, and a
// 								// pointer to it is stored within this object in the
// 								// corresponding place. The '1' indicates this image
// 								// will be used as a texture map. Use '2' to load
// 								// an image as a normal map, and '3' to load an
// 								// alpha map. Texture and normal maps are RGB .ppm
// 								// files, alpha maps are grayscale .pgm files.
// 								// * DO NOT * try to free image data loaded in this
// 								// way, the cleanup function already provided will do
// 								// this at the end.

//   insertObject(o,&object_list);			// <-- If you don't insert the object into the object list,
// 						//     nothing happens! your object won't be rendered.

// // normal mapping
// mirror at back
//  o=newPlane(.95,.05,.05,0.95,0.9,0.9,0.9,1,1,6);
// //   Scale(o,5,5,5);
// Translate(o,0,0,0.3);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newCyl(.05,.95,.95,.05,1,1,1,1,1,6);
//   RotateZ(o,-PI/2);
//   RotateX(o,PI/2);
//   RotateY(o,PI/2.5);
//   Scale(o,0.25,0.25,0.25);
// //   RotateZ(o,PI/3);
//  Translate(o,0,0,-0.5);
//  loadTexture(o,"./Texture/pockeball.ppm",1,&texture_list);  
//  loadTexture(o,"./Texture/npockeball.ppm",2,&texture_list);  
//   invert(&o->T[0][0],&o->Tinv[0][0]);
//   insertObject(o,&object_list);

//  o=newSphere(.05,.95,.95,0,1,.95,1,1,1,6);
//  Scale(o,0.25,0.25,0.25);
// //  Translate(o,1,0,1);
// //  RotateY(o,PI/2);
// //  RotateX(o,PI);
// //  RotateZ(o,-PI/2);
//  loadTexture(o,"./Texture/pockeball.ppm",1,&texture_list);  
// //  loadTexture(o,"./Texture/npockeball.ppm",2,&texture_list);  
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
//  Scale(o,0.25,0.25,0.25);
//  RotateX(o,PI/2);
//  RotateZ(o,PI/2);
//  Translate(o,-0.2, 0, -0.5);
//  loadTexture(o,"./Texture/pockeball.ppm",1,&texture_list);  
//  loadTexture(o,"./Texture/npockeball.ppm",2,&texture_list);  
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);



//  normal mapping LS
// o=newSphere(.05,.95,.95,.75,1,1,1,1,1,6);
// Scale(o,0.1,0.1,0.1);
// // Translate(o,-3,5,-3);
// Translate(o,0,5,0);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// o->isLightSource = 1;
// insertObject(o,&object_list);



// // refraction 1
//  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
//  Scale(o,5,5,5);
//  RotateX(o,PI/2);
//  Translate(o,0,-5,0);
//  loadTexture(o,"./Texture/ocean.ppm",1,&texture_list);  
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
//Scale(o,7, 7, 7);
//RotateZ(o,PI);
//Translate(o,0,2,6);
//loadTexture(o,"./Texture/luffy.ppm",1,&texture_list);
//invert(&o->T[0][0],&o->Tinv[0][0]);
//insertObject(o,&object_list);


//  o=newSphere(.05,.95,.95,.55,1,1,1,0,1.5,6);
//  Scale(o,.5,.5,.5);
//  Translate(o,0,0,-2);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);


// photon
  // o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
  // Scale(o,2,2,2);
  // RotateX(o,PI/2);
  // Translate(o,0,2,0);
  // loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
  // invert(&o->T[0][0],&o->Tinv[0][0]);
  // insertObject(o,&object_list);

  // o=newSphere(.05,.95,.95,.55,1,1,1,1,1.5,6);
  // Scale(o,2,2,2);
  // RotateY(o,-PI/2);
  // Translate(o,0,4,0);
  // loadTexture(o,"./Texture/pockeball.ppm",1,&texture_list);
  // loadTexture(o,"./Texture/npockeball.ppm",2,&texture_list);
  // invert(&o->T[0][0],&o->Tinv[0][0]);
  // insertObject(o,&object_list);



// Ligt source from hiiiiiiiigh
// o=newSphere(.05,.95,.95,.75,1,1,1,1,1,6);
// Scale(o,0.5,0.5,0.5);
// Translate(o,0,10,-10);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// o->isLightSource = 1;
// insertObject(o,&object_list);

// o=newPlane(.05,.95,.95,.75,1,1,1,1,1,6);
// RotateX(o,PI/2);
// Translate(o,0,5,1);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// o->isLightSource = 1;
// insertObject(o,&object_list);

// another angle for refraction
//  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
//  Scale(o,7,7,7);
//  Translate(o,7,0,9);
//  RotateZ(o,PI);
//  RotateY(o,PI/1.5);
//  loadTexture(o,"./Texture/luffy.ppm",1,&texture_list);  
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);


//  o=newSphere(.05,.95,.95,0.95,.75,.95,.55,0,1.5,6);
//  Scale(o,1.5,1.5,1);
//  Translate(o,0,-0.5,1);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

// Insert a single point light source. We set up its position as a point structure, and specify its
//  colour in terms of RGB (in [0,1]).
//  p.px=0;
//  p.py=25.5;
//  p.pz=-3.5;
//  p.pw=1;
//  l=newPLS(&p,.95,.95,.95);
//  insertPLS(l,&light_list);

// End of simple scene for Assignment 2
// Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
// or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
//
// Remember: A lot of the quality of your scene will depend on how much care you have put into defining
//           the relflectance properties of your objects, and the number and type of light sources
//           in the scene.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO: For Assignment 3 you *MUST* define your own cool scene.
//	   We will be looking for the quality of your scene setup, the use of hierarchical or composite
//	   objects that are more interesting than the simple primitives from A2, the use of textures
//        and other maps, illumination and illumination effects such as soft shadows, reflections and
//        transparency, and the overall visual quality of your result. Put some work into thinking
//        about these elements when designing your scene.
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
// Scale(o,2,2,2);
// RotateX(o,PI/2);
// Translate(o,0,-4,0);
// loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// insertObject(o,&object_list);

// o=newPlane(.05,.95,.95,.75,1,1,1,1,1,6);
// RotateX(o,PI/2);
// Translate(o,0,5,1);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// o->isLightSource = 1;
// insertObject(o,&object_list);

//ls
// o=newPlane(.05,.95,.95,.75,1,1,1,1,1,6);
// RotateX(o,PI/2);
// Translate(o,0,-7,0);
// invert(&o->T[0][0],&o->Tinv[0][0]);
// o->isLightSource = 1;
// insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.75,1,1,1,1,1,6);
Scale(o,0.25, 0.25, 0.25);
Translate(o,0,2,-1.9);
invert(&o->T[0][0],&o->Tinv[0][0]);
o->isLightSource = 1;
insertObject(o,&object_list);

//OBJ
//mirrors
o=newPlane(.95,0,0,0.6,1,1,1,1,1,6);
Scale(o,4, 4, 4);
RotateY(o,PI/6);
Translate(o,2*sqrt(3), 0, 2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newPlane(.95,0,0,0.6,1,1,1,1,1,6);
Scale(o,4, 4, 4);
RotateY(o,-PI/6);
Translate(o,-2*sqrt(3), 0, 2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

// background
o=newPlane(.05,.75,0,0,.55,.8,.75,1,1,2);
Scale(o,4, 4, 4);
RotateY(o,PI/6);
Translate(o,-2*sqrt(3), 0, -2);
loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newPlane(.05,.75,0,0,.55,.8,.75,1,1,2);
Scale(o,4, 4, 4);
RotateY(o,-PI/6);
Translate(o,2*sqrt(3), 0, -2);
loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

// bottom
o=newPlane(.05,.75,0,0,.55,.8,.75,1,1,2);
RotateX(o,PI/2);
Translate(o,0,-4,0);
RotateY(o,PI/4);
Scale(o,2*sqrt(6),1,2*sqrt(2));
loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

// top
o=newPlane(.05,.75,0,0,.55,.8,.75,1,1,2);
RotateX(o,PI/2);
Translate(o,0,4,0);
RotateY(o,PI/4);
Scale(o,2*sqrt(6),1,2*sqrt(2));
loadTexture(o,"./Texture/universe.ppm",1,&texture_list);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

// photon ball
 o=newSphere(.05,.95,.95,.95,1,1,1,0,1.5,6);
 Scale(o,0.75, 0.75, 0.75);
 Translate(o,0,-3,0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // photon frame
 o=newPlane(.05,0.75,0,0,.55,.8,.75,1,1,2);
 RotateZ(o,PI);
 Scale(o,0.85,1,1); 
 RotateY(o,PI);
 Translate(o,0, -1, -3.29);
 loadTexture(o,"./Texture/paco.ppm",1,&texture_list);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(.05,.75,0,0,.55,.8,.75,1,1,2);
 Scale(o,0.85*1.3,1*1.3,1);
 RotateY(o,PI);
 Translate(o,0, -1, -3.3);
 loadTexture(o,"./Texture/frame.ppm",1,&texture_list);
 loadTexture(o,"./Texture/nframe.ppm",2,&texture_list);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //planet
 o=newSphere(.05,.95,.95,.55,1,1,1,1,1.5,6);
 Scale(o,0.5, 0.5, 0.5);
 RotateX(o,PI/2);
 Translate(o,-1,1,-1.5);
 loadTexture(o,"./Texture/planet1.PPM",1,&texture_list);
 loadTexture(o,"./Texture/nplanet1.PPM",2,&texture_list);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.95,.55,1,1,1,1,1.5,6);
 Scale(o,0.75, 0.75, 0.75);
 RotateX(o,PI/4);
 Translate(o,1.5,1.5,-1);
 loadTexture(o,"./Texture/planet2.PPM",1,&texture_list);
 loadTexture(o,"./Texture/nplanet2.PPM",2,&texture_list);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);