/*
   CSC D18 - Path Tracer
   
   Utilities for the Path Tracer. 
   
   Derived from the Ray Tracer code. You should already be familiar 
   with most of what is here.
   
   Last updated: Aug. 2017  -  F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:  Zhou Xiaohan 
* 2) Student Name:  Hu Yuxuan
*
* 1) Student number:  1000763836
* 2) Student number:  1000351460 
*
* 1) UtorID  zhouxi47
* 2) UtorID  huyuxuan
*
* We hereby certify that the work contained here is our own
*
* ______ZhouXiaoHan_____             ______huyuxuanKevin_________
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils_path.h"

// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4]={{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point3D *newPoint(double px, double py, double pz)
{
 // Allocate a new point structure, initialize it to
 // the specified coordinates, and return a pointer
 // to it.

 struct point3D *pt=(struct point3D *)calloc(1,sizeof(struct point3D));
 if (!pt) fprintf(stderr,"Out of memory allocating point structure!\n");
 else
 {
  pt->px=px;
  pt->py=py;
  pt->pz=pz;
  pt->pw=1.0;
 }
 return(pt);
}


/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{
 // Transforms a ray using the inverse transform for the specified object. This is so that we can
 // use the intersection test for the canonical object. Note that this has to be done carefully!

 ///////////////////////////////////////////
 // TO DO: Complete this function
 ///////////////////////////////////////////
 memcpy(ray_transformed, ray_orig, sizeof(struct ray3D));
 
 // transfer ray->p0 to canonical object
     matVecMult(obj->Tinv, &ray_transformed->p0);
 
 // transfer ray->d to canonical object
     ray_transformed->d.pw = 0;
     matVecMult(obj->Tinv, &ray_transformed->d);
 
     ray_transformed->d.pw = 1;
}

inline void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
 // Computes the normal at an affinely transformed point given the original normal and the
 // object's inverse transformation. From the notes:
 // n_transformed=A^-T*n normalized.

 ///////////////////////////////////////////
 // TO DO: Complete this function
 ///////////////////////////////////////////
 memcpy(n_transformed, n_orig, sizeof(struct point3D));
 
     double At[4][4];
     for (int i = 0; i < 4; i++) {
         for (int j = 0; j < 4; j++) {
             At[i][j] = obj->Tinv[j][i];
         }
     }
 
 // transfer from canonical object back to the specific object
     n_transformed->pw = 0;
     matVecMult(At, n_transformed);
     normalize(n_transformed);
     n_transformed->pw = 1;
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
void insertObject(struct object3D *o, struct object3D **list)
{
 if (o==NULL) return;
 // Inserts an object into the object list.
 if (*(list)==NULL)
 {
  *(list)=o;
  (*(list))->next=NULL;
 }
 else
 {
  o->next=(*(list))->next;
  (*(list))->next=o;
 }
}

struct object3D *newPlane(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
 // Intialize a new plane with the specified parameters:
 // diffPct, reflPct, tranPct - specify the amount of diffuse, reflective, and
 //   refracting properties of the material. They *must* sum to 1.0 
 // r, g, b, - Colour for this plane
 // refl_sig - Determines the amount of spread for reflection directions. If zero
 //   rays are reflected only along the perfect reflection direction, if non-zero,
 //   the perfect reflection direction is bent a bit (the amount is drawn from a
 //   zero-mean Gaussian distribution with sigma refl_sig). This makes the reflection
 //   component less sharp, and makes the material look more 'matte'
 // r_index - Refraction index for the refraction component.
 //
 // The plane is defined by the following vertices (CCW)
 // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
 // With normal vector (0,0,1) (i.e. parallel to the XY plane)

 struct object3D *plane=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!plane) fprintf(stderr,"Unable to allocate new plane, out of memory!\n");
 else
 {
  plane->diffPct=diffPct;
  plane->reflPct=reflPct;
  plane->tranPct=tranPct;
  plane->col.R=r;
  plane->col.G=g;
  plane->col.B=b;
  plane->refl_sig=refl_sig;
  plane->r_index=r_index;
  plane->intersect=&planeIntersect;
  plane->surfaceCoords=&planeCoordinates;
  plane->randomPoint=&planeSample;
  plane->texImg=NULL;
  plane->normalMap=NULL;
  memcpy(&plane->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&plane->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  plane->textureMap=&texMap;
  plane->frontAndBack=1;
  plane->normalMapped=0;
  plane->isCSG=0;
  plane->isLightSource=0;
  plane->LSweight=1.0;
  plane->CSGnext=NULL;
  plane->next=NULL;
}
 return(plane);
}

struct object3D *newSphere(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
 // Intialize a new sphere with the specified parameters. The parameters have the same meaning
 // as for planes, so have a look above at the newPlane() function comments.
 //
 // This is assumed to represent a unit sphere centered at the origin.

 struct object3D *sphere=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!sphere) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  sphere->diffPct=diffPct;
  sphere->reflPct=reflPct;
  sphere->tranPct=tranPct;
  sphere->col.R=r;
  sphere->col.G=g;
  sphere->col.B=b;
  sphere->refl_sig=refl_sig;
  sphere->r_index=r_index;
  sphere->intersect=&sphereIntersect;
  sphere->surfaceCoords=&sphereCoordinates;
  sphere->randomPoint=&sphereSample;
  sphere->texImg=NULL;
  sphere->normalMap=NULL;
  memcpy(&sphere->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&sphere->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  sphere->textureMap=&texMap;
  sphere->frontAndBack=0;
  sphere->normalMapped=0;
  sphere->isCSG=0;
  sphere->isLightSource=0;
  sphere->LSweight=1.0;
  sphere->CSGnext=NULL;
  sphere->next=NULL; }
 return(sphere);
}

struct object3D *newCyl(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
 ///////////////////////////////////////////////////////////////////////////////////////
 // TO DO:
 //	Complete the code to create and initialize a new cylinder object.
 ///////////////////////////////////////////////////////////////////////////////////////  
  
 struct object3D *cylinder=(struct object3D *)calloc(1,sizeof(struct object3D));
 
  if (!cylinder) fprintf(stderr,"Unable to allocate new cylinder, out of memory!\n");
  else
  {
    cylinder->diffPct=diffPct;
    cylinder->reflPct=reflPct;
    cylinder->tranPct=tranPct;
    cylinder->col.R=r;
    cylinder->col.G=g;
    cylinder->col.B=b;
    cylinder->refl_sig=refl_sig;
    cylinder->r_index=r_index;
    cylinder->intersect=&cylIntersect;
    cylinder->surfaceCoords=&cylCoordinates;
    cylinder->randomPoint=&cylSample;
    cylinder->texImg=NULL;
    cylinder->normalMap=NULL;
    memcpy(&cylinder->T[0][0],&eye4x4[0][0],16*sizeof(double));
    memcpy(&cylinder->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
    cylinder->textureMap=&texMap;
    cylinder->frontAndBack=0;
    cylinder->normalMapped=0;
    cylinder->isCSG=0;
    cylinder->isLightSource=0;
    cylinder->LSweight=1.0;
    cylinder->CSGnext=NULL;
    cylinder->next=NULL; }
    return(cylinder);
}


///////////////////////////////////////////////////////////////////////////////////////
// TO DO:
//	Complete the functions that compute intersections for the canonical plane
//      and canonical sphere with a given ray. This is the most fundamental component
//      of the raytracer.
///////////////////////////////////////////////////////////////////////////////////////
void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical plane.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
// transfer ray to canonical coordinates
struct ray3D ray_transformed;
rayTransform(ray, &ray_transformed, plane);

struct point3D plane_n;
plane_n.px = 0;
plane_n.py = 0;
plane_n.pz = 1;

double temp = dot(&ray_transformed.d, &plane_n);
if (temp > 0) {
    // angle between ray->d and normal < 90
    plane_n.px *= -1;
    plane_n.py *= -1;
    plane_n.pz *= -1;
    temp = dot(&ray_transformed.d, &plane_n);
}

struct point3D p1;
p1.px = -ray_transformed.p0.px;
p1.py = -ray_transformed.p0.py;
p1.pz = -ray_transformed.p0.pz;

*lambda = dot(&p1, &plane_n) / temp;
(ray_transformed.rayPos)(&ray_transformed, *lambda, p);

if (fabs(p->px) > 1 || fabs(p->py) > 1) {
    // point not on the plane
    *lambda = -1;
} else {
    // a and b
    *a = (p->px + 1) / 2;
    *b = (p->py + 1) / 2;

    // bump mapping
    if (plane->normalMap != NULL) {
        double bump_R, bump_G, bump_B;
        plane->textureMap(plane->normalMap, *a, *b, &bump_R, &bump_G, &bump_B);

        plane_n.px = 2 * bump_R - 1;
        plane_n.py = 2 * bump_B - 1;
        plane_n.pz = 2 * bump_G - 1;

        normalize(&plane_n);
    }

    // transfer intersection point to specific coordinate
    (ray->rayPos)(ray, *lambda, p);

    // transfer plane's normal to specific coordinate
    normalTransform(&plane_n, n, plane);
}

}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical sphere.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
// transfer ray to canonical coordinates
struct ray3D ray_transformed;
rayTransform(ray, &ray_transformed, sphere);

// point a is ray.po, point d is ray.d
double cap_a = dot(&ray_transformed.d, &ray_transformed.d);
double cap_b = dot(&ray_transformed.p0, &ray_transformed.d);
double cap_c = dot(&ray_transformed.p0, &ray_transformed.p0) - 1;
double cap_d = cap_b * cap_b - cap_a * cap_c;

if (cap_d < 0) // no interesection
{
    *lambda = -1;
} else {
    double lambda1 = -cap_b / cap_a + sqrt(cap_d) / cap_a;
    double lambda2 = -cap_b / cap_a - sqrt(cap_d) / cap_a;

    // cap_d == 0, one intersction, ray grazes the sphere
    // l1>l2 && l2>0, both in front of plane, l2 is the closest.
    *lambda = lambda2;
    if (lambda1 < 0 && lambda2 < 0) // both behind the view-plane, not visible
    {
        *lambda = -1;
    } else if (lambda1 > 0 && lambda2 < 0) // only l1 is visile
    {
        *lambda = lambda1;
    }
}

if (*lambda >= 0) {
    (ray_transformed.rayPos)(&ray_transformed, *lambda, p);
    struct point3D sphere_n;
    sphere_n.px = p->px;
    sphere_n.py = p->py;
    sphere_n.pz = p->pz;

    if (dot(&sphere_n, &ray_transformed.d) > 0) {
        sphere_n.px *= -1;
        sphere_n.py *= -1;
        sphere_n.pz *= -1;
    }

    // a and b
    double phi = p->px != 0 ? atan(p->py/p->px) : PI/2;
    if (p->px < 0) {
        phi += PI;
    } else if (p->py < 0) {
        phi += 2*PI;
    }
    *a = phi/(2*PI);

    *b = (-asin(p->pz) + PI / 2) / PI;

    // bump mapping
    if (sphere->normalMap) {
        double bump_R, bump_G, bump_B;
        sphere->textureMap(sphere->normalMap, *a, *b, &bump_R, &bump_G, &bump_B);

        struct point3D normalTex;
        normalTex.px = 2 * bump_R - 1;
        normalTex.py = 2 * bump_G - 1;
        normalTex.pz = 2 * bump_B - 1;
        struct point3D tangent;
        tangent.px = sphere_n.py;
        tangent.py = -sphere_n.px;
        tangent.pz = 0;
        struct point3D *binormal = cross(&tangent, &sphere_n);

        double M[4][4] = {
                {sphere_n.py,  binormal->px, sphere_n.px, 0.0},
                {-sphere_n.px, binormal->py, sphere_n.py, 0.0},
                {sphere_n.pz,  binormal->pz, sphere_n.pz, 0.0},
                {0.0,           0.0,          0.0,          1.0}};

        matVecMult(M, &normalTex);
        normalize(&normalTex);

        sphere_n.px = normalTex.px;
        sphere_n.py = normalTex.py;
        sphere_n.pz = normalTex.pz;

        free(binormal);
    }

    // transfer intersection point to specific coordinate
    (ray->rayPos)(ray, *lambda, p);

    // transfer plane's normal to specific coordinate
    normalTransform(&sphere_n, n, sphere);

}
}

void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Computes and returns the value of 'lambda' at the intersection
 // between the specified ray and the specified canonical cylinder.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 struct ray3D ray_transformed;
 rayTransform(r, &ray_transformed, cylinder);

 double cap_a = dot(&ray_transformed.d, &ray_transformed.d) - pow(ray_transformed.d.pz, 2);
 double cap_b = dot(&ray_transformed.p0, &ray_transformed.d) - ray_transformed.p0.pz * ray_transformed.d.pz;
 double cap_c = dot(&ray_transformed.p0, &ray_transformed.p0) - pow(ray_transformed.p0.pz, 2) - 1;
 double cap_d = cap_b * cap_b - cap_a * cap_c;

 if (cap_d < 0) // no interesection
 {
     *lambda = -1;
 } else {
     double lambda1 = -cap_b / cap_a + sqrt(cap_d) / cap_a;
     double lambda2 = -cap_b / cap_a - sqrt(cap_d) / cap_a;

     // cap_d == 0, one intersction, ray grazes the sphere
     // l1>l2 && l2>0, both in front of plane, l2 is the closest.
     *lambda = lambda2;
     (ray_transformed.rayPos)(&ray_transformed, *lambda, p);
     if ((lambda1 > 0 && lambda2 < 0) || (p->pz > 0.5 || p->pz < -0.5)) {
         // only l1 is visile || pz out of range
         *lambda = lambda1;
     }

     (ray_transformed.rayPos)(&ray_transformed, *lambda, p);
     if ((lambda1 < 0 && lambda2 < 0) || (p->pz > 0.5 || p->pz < -0.5)) {
         // both behind the view-plane, not visible
         // if pz for both l1 and l2 out of range
         *lambda = -1;
     }
 }

 if (*lambda >= 0) {
     (ray_transformed.rayPos)(&ray_transformed, *lambda, p);

     struct point3D cylinder_n;
     cylinder_n.px = p->px;
     cylinder_n.py = p->py;
     cylinder_n.pz = 0;
     if (dot(&cylinder_n, &ray_transformed.d) > 0) {
         cylinder_n.px *= -1;
         cylinder_n.py *= -1;
         cylinder_n.pz *= -1;
     }

     // a and b
    //  *a = (acos(p->px) + PI / 2) / PI;
    //  *b = p->pz + 0.5;
    double phi = p->px != 0 ? atan(p->py/p->px) : PI/2;
    if (p->px < 0) {
        phi += PI;
    } else if (p->py < 0) {
        phi += 2*PI;
    }
    *a = phi/(2*PI);
  *b = p->pz + 0.5;

// bump mapping
     if (cylinder->normalMap) {
         double bump_R, bump_G, bump_B;
         cylinder->textureMap(cylinder->normalMap, *a, *b, &bump_R, &bump_G, &bump_B);

         struct point3D *normalTex = newPoint(2 * bump_R - 1, 2 * bump_G - 1, 2 * bump_B - 1);
         struct point3D *tangent = newPoint(cylinder_n.py, -cylinder_n.px, 0);
         struct point3D *binormal = cross(tangent, &cylinder_n);

         double M[4][4] = {
                 {cylinder_n.py,  binormal->px, cylinder_n.px, 0.0},
                 {-cylinder_n.px, binormal->py, cylinder_n.py, 0.0},
                 {cylinder_n.pz,  binormal->pz, cylinder_n.pz, 0.0},
                 {0.0,             0.0,          0.0,            1.0}};

         matVecMult(M, normalTex);
         normalize(normalTex);

         cylinder_n.px = normalTex->px;
         cylinder_n.py = normalTex->py;
         cylinder_n.pz = normalTex->pz;

         free(normalTex);
         free(binormal);
         free(tangent);
     }

     // transfer intersection point to specific coordinate
     (r->rayPos)(r, *lambda, p);

     // transfer cylinder's normal to specific coordinate
     normalTransform(&cylinder_n, n, cylinder);
 }  
}

/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z)
{
 // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
 // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
 // bottom of the plane.
 
 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 struct point3D *p = newPoint(2 * a - 1, 2 * b - 1, 0);
 matVecMult(plane->T, p);
 *x = p->px;
 *y = p->py;
 *z = p->pz;
 free(p);   
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x, double *y, double *z)
{
 // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
 // 'a' in [0, 2*PI] corresponds to the spherical coordinate theta
 // 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////   
 double u = cos(b);
 double temp = sqrt(1 - pow(u, 2));

 struct point3D *p = newPoint(temp * cos(a), temp * sin(a), u);
 matVecMult(sphere->T, p);
 *x = p->px;
 *y = p->py;
 *z = p->pz;
 free(p);
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x, double *y, double *z)
{
 // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
 // 'a' in [0, 2*PI] corresponds to angle theta around the cylinder
 // 'b' in [0, 1] corresponds to height from the bottom
 
 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////  
 struct point3D *p = newPoint(cos(a), sin(a), b);
 
     matVecMult(cyl->T, p);
     *x = p->px;
     *y = p->py;
     *z = p->pz;
     free(p); 
}

void planeSample(struct object3D *plane, double *x, double *y, double *z)
{
 // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
 // Sapling should be uniform, meaning there should be an equal change of gedtting
 // any spot on the plane

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////   
 double a = (double) rand() / (double) RAND_MAX;
 double b = (double) rand() / (double) RAND_MAX;
 (plane->surfaceCoords)(plane, a, b, x, y, z);
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z)
{
 // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the sphere
 // Sampling should be uniform - note that this is tricky for a sphere, do some
 // research and document in your report what method is used to do this, along
 // with a reference to your source.
 
 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 // references http://mathworld.wolfram.com/SpherePointPicking.html
// 'a' in [0, 2PI], 'b' in [-PI/2, PI/2]
 double a = 2 * PI * (double) rand() / RAND_MAX;
 double b = (PI * (double) rand() / RAND_MAX) - PI / 2;

 (sphere->surfaceCoords)(sphere, a, b, x, y, z);  
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z)
{
 // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the cylinder
 // Sampling should be uniform over the cylinder.

 /////////////////////////////////
 // TO DO: Complete this function.
 /////////////////////////////////
 double b = 2 * PI * (double) rand() / RAND_MAX;
 double a = ((double) rand() / RAND_MAX) - 0.5;

 (cyl->surfaceCoords)(cyl, a, b, x, y, z);   
}

//////////////////////////////////
// Importance sampling for BRDF
//////////////////////////////////
void cosWeightedSample(struct point3D *n, struct point3D *d)
{
 // This function returns a randomly sampled direction over
 // a hemisphere whose pole is the normal direction n. The
 // sampled direction comes from a distribution weighted
 // by the cosine of the angle between n and d.
 // Use this for importance sampling for diffuse surfaces.
  
 double u1,r,theta,phi;
 double x,y,z,c;
 double v[4][4],R[4][4];
 struct point3D nz,*cr;
 char line[1024];
 
 // Random sample on hemisphere with cosine-weighted distribution
 u1=drand48();
 r=sqrt(u1);
 theta=2*PI*drand48();
 x=r*cos(theta);
 y=r*sin(theta);
 z=sqrt(1.0-(x*x)-(y*y));
 
 // Need a rotation matrix - start with identity
 memset(&R[0][0],0,4*4*sizeof(double));
 R[0][0]=1.0;
 R[1][1]=1.0;
 R[2][2]=1.0;
 R[3][3]=1.0;

 // Rotation based on cylindrical coordinate conversion
 theta=atan2(n->py,n->px);
 phi=acos(n->pz);
 RotateYMat(R,phi);
 RotateZMat(R,theta);

 // Rotate d to align with normal 
 d->px=x;
 d->py=y;
 d->pz=z;
 d->pw=1.0;
 matVecMult(R,d);
 
 return;
}

/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list)
{
 // Load a texture or normal map image from file and assign it to the
 // specified object. 
 // type:   1  ->  Texture map  (RGB, .ppm)
 //         2  ->  Normal map   (RGB, .ppm)
 //         3  ->  Alpha map    (grayscale, .pgm)
 // Stores loaded images in a linked list to avoid replication
 struct image *im;
 struct textureNode *p;
 
 if (o!=NULL)
 {
  // Check current linked list
  p=*(t_list);
  while (p!=NULL)
  {
   if (strcmp(&p->name[0],filename)==0)
   {
    // Found image already on the list
    if (type==1) o->texImg=p->im;
    else if (type==2) o->normalMap=p->im;
    else o->alphaMap=p->im;
    return;
   }
   p=p->next;
  }    

  // Load this texture image 
  if (type==1||type==2)
   im=readPPMimage(filename);
  else if (type==3)
   im=readPGMimage(filename);

  // Insert it into the texture list
  if (im!=NULL)
  {
   p=(struct textureNode *)calloc(1,sizeof(struct textureNode));
   strcpy(&p->name[0],filename);
   p->type=type;
   p->im=im;
   p->next=NULL; 
   // Insert into linked list
   if ((*(t_list))==NULL)
    *(t_list)=p;
   else
   {
    p->next=(*(t_list))->next;
    (*(t_list))->next=p;
   }
   // Assign to object
   if (type==1) o->texImg=im;
   else if (type==2) o->normalMap=im;
   else o->alphaMap=im;
  }
 
 }  // end if (o != NULL)
}


void texMap(struct image *img, double a, double b, double *R, double *G, double *B)
{
 /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */

 //////////////////////////////////////////////////
 // TO DO:
 //
 //  Complete this function to return the colour
 // of the texture image at the specified texture
 // coordinates. Your code should use bi-linear
 // interpolation to obtain the texture colour.
 //////////////////////////////////////////////////
 int ia, ib, ia2, ib2;     // Coordinates of interpolation vertices
 double r1, g1, b1;       // Colour at the four vertex pixels
 double r2, g2, b2;
 double r3, g3, b3;
 double r4, g4, b4;
 double fracc_a, fracc_b;               // Location within 4-pixel region
 double *ip = (double *) img->rgbdata;     // Pointer to double image data

// Obtain integer coordinates in terms of texture size. Ensure they will
// be within bounds
 if (a < 0) a = 0;
 if (b < 0) b = 0;
 if (a > 1) a = 1;
 if (b > 1) b = 1;

 fracc_a = (a * img->sx);
 fracc_b = (b * img->sy);
 ia = (int) fracc_a;
 ib = (int) fracc_b;
 if (ia >= img->sx) ia = img->sx - 1;
 if (ib >= img->sy) ib = img->sy - 1;

 fracc_a = fracc_a - ia;
 fracc_b = fracc_b - ib;
 if (fracc_a > 0) ia2 = ia + 1; else ia2 = ia;
 if (fracc_b > 0) ib2 = ib + 1; else ib2 = ib;

 if (ia2 >= img->sx) ia2 = img->sx - 1;
 if (ib2 >= img->sy) ib2 = img->sy - 1;

 r1 = *(ip + ((ia + (ib * img->sx)) * 3) + 0);
 g1 = *(ip + ((ia + (ib * img->sx)) * 3) + 1);
 b1 = *(ip + ((ia + (ib * img->sx)) * 3) + 2);

 r2 = *(ip + ((ia2 + (ib * img->sx)) * 3) + 0);
 g2 = *(ip + ((ia2 + (ib * img->sx)) * 3) + 1);
 b2 = *(ip + ((ia2 + (ib * img->sx)) * 3) + 2);

 r3 = *(ip + ((ia + (ib2 * img->sx)) * 3) + 0);
 g3 = *(ip + ((ia + (ib2 * img->sx)) * 3) + 1);
 b3 = *(ip + ((ia + (ib2 * img->sx)) * 3) + 2);

 r4 = *(ip + ((ia2 + (ib2 * img->sx)) * 3) + 0);
 g4 = *(ip + ((ia2 + (ib2 * img->sx)) * 3) + 1);
 b4 = *(ip + ((ia2 + (ib2 * img->sx)) * 3) + 2);

 *R = (1 - fracc_b) * (((1 - fracc_a) * r1) + (fracc_a * r2));
 *R += fracc_b * (((1 - fracc_a) * r3) + (fracc_a * r4));

 *G = (1 - fracc_b) * (((1 - fracc_a) * g1) + (fracc_a * g2));
 *G += fracc_b * (((1 - fracc_a) * g3) + (fracc_a * g4));

 *B = (1 - fracc_b) * (((1 - fracc_a) * b1) + (fracc_a * b2));
 *B += fracc_b * (((1 - fracc_a) * b3) + (fracc_a * b4));

 return;
}

void alphaMap(struct image *img, double a, double b, double *alpha)
{
 // Just like texture map but returns the alpha value at a,b,
 // notice that alpha maps are single layer grayscale images, hence
 // the separate function.

 //////////////////////////////////////////////////
 // TO DO:
 //
 //  Complete this function to return the alpha
 // value from the image at the specified texture
 // coordinates. Your code should use bi-linear
 // interpolation to obtain the texture colour.
 //////////////////////////////////////////////////
 
 if (a<0) a=0; if (b<0) b=0; if (a>1) a=1; if (b>1) b=1;
 double newa, newb;
 int ia, ib;
 double *imgcol = (double *)img->rgbdata;
 
 newa = a * (img->sx-1);
 newb = b * (img->sy-1);
 ia = (int) newa;
 ib = (int) newb;
 
 //printf("im here\n");
 *(alpha) = *(imgcol + ia + (ib*img->sx)); 
}


///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv)
{
 // Computes the inverse of transformation matrix T.
 // the result is returned in Tinv.

 double *U, *s, *V, *rv1;
 int singFlag, i;

 // Invert the affine transform
 U=NULL;
 s=NULL;
 V=NULL;
 rv1=NULL;
 singFlag=0;

 SVD(T,4,4,&U,&s,&V,&rv1);
 if (U==NULL||s==NULL||V==NULL)
 {
  fprintf(stderr,"Error: Matrix not invertible for this object, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Check for singular matrices...
 for (i=0;i<4;i++) if (*(s+i)<1e-9) singFlag=1;
 if (singFlag)
 {
  fprintf(stderr,"Error: Transformation matrix is singular, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 // Compute and store inverse matrix
 InvertMatrix(U,s,V,4,Tinv);

 free(U);
 free(s);
 free(V);
}

void RotateXMat(double T[4][4], double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // X axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateX(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // X axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateYMat(double T[4][4], double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Y axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateY(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Y axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateZMat(double T[4][4], double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Z axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateZ(struct object3D *o, double theta)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that rotates the object theta *RADIANS* around the
 // Z axis.

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,o->T);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that translates the object by the specified amounts.

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that translates the object by the specified amounts.

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that scales the object as indicated.

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{
 // Multiply the current object transformation matrix T in object o
 // by a matrix that scales the object as indicated.

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,o->T);
 o->LSweight*=(sx*sy*sz);	// Update object volume! careful
				// won't work for hierarchical
				// objects!
}

void printmatrix(double mat[4][4])
{
 fprintf(stderr,"Matrix contains:\n");
 fprintf(stderr,"%f %f %f %f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize)
{
 /*
   This function sets up the camera axes and viewing direction as discussed in the
   lecture notes.
   e - Camera center
   g - Gaze direction
   up - Up vector
   fov - Fild of view in degrees
   f - focal length
 */
 struct view *c;
 struct point3D *u, *v;

 u=v=NULL;

 // Allocate space for the camera structure
 c=(struct view *)calloc(1,sizeof(struct view));
 if (c==NULL)
 {
  fprintf(stderr,"Out of memory setting up camera model!\n");
  return(NULL);
 }

 // Set up camera center and axes
 c->e.px=e->px;		// Copy camera center location, note we must make sure
 c->e.py=e->py;		// the camera center provided to this function has pw=1
 c->e.pz=e->pz;
 c->e.pw=1;

 // Set up w vector (camera's Z axis). w=-g/||g||
 c->w.px=-g->px;
 c->w.py=-g->py;
 c->w.pz=-g->pz;
 c->w.pw=1;
 normalize(&c->w);

 // Set up the horizontal direction, which must be perpenticular to w and up
 u=cross(&c->w, up);
 normalize(u);
 c->u.px=u->px;
 c->u.py=u->py;
 c->u.pz=u->pz;
 c->u.pw=1;

 // Set up the remaining direction, v=(u x w)  - Mind the signs
 v=cross(&c->u, &c->w);
 normalize(v);
 c->v.px=v->px;
 c->v.py=v->py;
 c->v.pz=v->pz;
 c->v.pw=1;

 // Copy focal length and window size parameters
 c->f=f;
 c->wl=wl;
 c->wt=wt;
 c->wsize=wsize;

 // Set up coordinate conversion matrices
 // Camera2World matrix (M_cw in the notes)
 // Mind the indexing convention [row][col]
 c->C2W[0][0]=c->u.px;
 c->C2W[1][0]=c->u.py;
 c->C2W[2][0]=c->u.pz;
 c->C2W[3][0]=0;

 c->C2W[0][1]=c->v.px;
 c->C2W[1][1]=c->v.py;
 c->C2W[2][1]=c->v.pz;
 c->C2W[3][1]=0;

 c->C2W[0][2]=c->w.px;
 c->C2W[1][2]=c->w.py;
 c->C2W[2][2]=c->w.pz;
 c->C2W[3][2]=0;

 c->C2W[0][3]=c->e.px;
 c->C2W[1][3]=c->e.py;
 c->C2W[2][3]=c->e.pz;
 c->C2W[3][3]=1;

 // World2Camera matrix (M_wc in the notes)
 // Mind the indexing convention [row][col]
 c->W2C[0][0]=c->u.px;
 c->W2C[1][0]=c->v.px;
 c->W2C[2][0]=c->w.px;
 c->W2C[3][0]=0;

 c->W2C[0][1]=c->u.py;
 c->W2C[1][1]=c->v.py;
 c->W2C[2][1]=c->w.py;
 c->W2C[3][1]=0;

 c->W2C[0][2]=c->u.pz;
 c->W2C[1][2]=c->v.pz;
 c->W2C[2][2]=c->w.pz;
 c->W2C[3][2]=0;

 c->W2C[0][3]=-dot(&c->u,&c->e);
 c->W2C[1][3]=-dot(&c->v,&c->e);
 c->W2C[2][3]=-dot(&c->w,&c->e);
 c->W2C[3][3]=1;

 free(u);
 free(v);
 return(c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename)
{
 // Reads an image from a .ppm file. A .ppm file is a very simple image representation
 // format with a text header followed by the binary RGB data at 24bits per pixel.
 // The header has the following form:
 //
 // P6
 // # One or more comment lines preceded by '#'
 // 340 200
 // 255
 //
 // The first line 'P6' is the .ppm format identifier, this is followed by one or more
 // lines with comments, typically used to inidicate which program generated the
 // .ppm file.
 // After the comments, a line with two integer values specifies the image resolution
 // as number of pixels in x and number of pixels in y.
 // The final line of the header stores the maximum value for pixels in the image,
 // usually 255.
 // After this last header line, binary data stores the RGB values for each pixel
 // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
 //
 // NOTE: Windows file handling is rather crotchetty. You may have to change the
 //       way this file is accessed if the images are being corrupted on read
 //       on Windows.
 //
 // readPPMdata converts the image colour information to floating point. This is so that
 // the texture mapping function doesn't have to do the conversion every time
 // it is asked to return the colour at a specific location.
 //

 FILE *f;
 struct image *im;
 char line[1024];
 int sizx,sizy;
 int i;
 unsigned char *tmp;
 double *fRGB;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  f=fopen(filename,"rb+");
  if (f==NULL)
  {
   fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
   free(im);
   return(NULL);
  }
  fgets(&line[0],1000,f);
  if (strcmp(&line[0],"P6\n")!=0)
  {
   fprintf(stderr,"Wrong file format, not a .ppm file or header end-of-line characters missing\n");
   free(im);
   fclose(f);
   return(NULL);
  }
  fprintf(stderr,"%s\n",line);
  // Skip over comments
  fgets(&line[0],511,f);
  while (line[0]=='#')
  {
   fprintf(stderr,"%s",line);
   fgets(&line[0],511,f);
  }
  sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
  fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
  im->sx=sizx;
  im->sy=sizy;

  fgets(&line[0],9,f);  	                // Read the remaining header line
  fprintf(stderr,"%s\n",line);
  tmp=(unsigned char *)calloc(sizx*sizy*3,sizeof(unsigned char));
  fRGB=(double *)calloc(sizx*sizy*3,sizeof(double));
  if (tmp==NULL||fRGB==NULL)
  {
   fprintf(stderr,"Out of memory allocating space for image\n");
   free(im);
   fclose(f);
   return(NULL);
  }

  fread(tmp,sizx*sizy*3*sizeof(unsigned char),1,f);
  fclose(f);

  // Conversion to floating point
  for (i=0; i<sizx*sizy*3; i++) *(fRGB+i)=((double)*(tmp+i))/255.0;
  free(tmp);
  im->rgbdata=(void *)fRGB;

  return(im);
 }

 fprintf(stderr,"Unable to allocate memory for image structure\n");
 return(NULL);
}

struct image *readPGMimage(const char *filename)
{
 // Just like readPPMimage() except it is used to load grayscale alpha maps. In
 // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0 
 // correspondst to alpha=0 (fully transparent).
 // A .pgm header of the following form is expected:
 //
 // P5
 // # One or more comment lines preceded by '#'
 // 340 200
 // 255
 //
 // readPGMdata converts the image grayscale data to double floating point in [0,1]. 

 FILE *f;
 struct image *im;
 char line[1024];
 int sizx,sizy;
 int i;
 unsigned char *tmp;
 double *fRGB;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  f=fopen(filename,"rb+");
  if (f==NULL)
  {
   fprintf(stderr,"Unable to open file %s for reading, please check name and path\n",filename);
   free(im);
   return(NULL);
  }
  fgets(&line[0],1000,f);
  if (strcmp(&line[0],"P5\n")!=0)
  {
   fprintf(stderr,"Wrong file format, not a .pgm file or header end-of-line characters missing\n");
   free(im);
   fclose(f);
   return(NULL);
  }
  // Skip over comments
  fgets(&line[0],511,f);
  while (line[0]=='#')
   fgets(&line[0],511,f);
  sscanf(&line[0],"%d %d\n",&sizx,&sizy);           // Read file size
  im->sx=sizx;
  im->sy=sizy;

  fgets(&line[0],9,f);  	                // Read the remaining header line
  tmp=(unsigned char *)calloc(sizx*sizy,sizeof(unsigned char));
  fRGB=(double *)calloc(sizx*sizy,sizeof(double));
  if (tmp==NULL||fRGB==NULL)
  {
   fprintf(stderr,"Out of memory allocating space for image\n");
   free(im);
   fclose(f);
   return(NULL);
  }

  fread(tmp,sizx*sizy*sizeof(unsigned char),1,f);
  fclose(f);

  // Conversion to double floating point
  for (i=0; i<sizx*sizy; i++) *(fRGB+i)=((double)*(tmp+i))/255.0;
  free(tmp);
  im->rgbdata=(void *)fRGB;

  return(im);
 }

 fprintf(stderr,"Unable to allocate memory for image structure\n");
 return(NULL);
}

struct image *newImage(int size_x, int size_y)
{
 // Allocates and returns a new image with all zeros. This allocates a double
 // precision floating point image! MIND the difference with the raytracer
 // code that uses 24bpp images.
 struct image *im;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  im->sx=size_x;
  im->sy=size_y;
  im->rgbdata=(void *)calloc(size_x*size_y*3,sizeof(double));
  if (im->rgbdata!=NULL) return(im);
 }
 fprintf(stderr,"Unable to allocate memory for new image\n");
 return(NULL);
}

void imageOutput(struct image *im, const char *filename)
{
 // Writes out a .ppm file from the image data contained in 'im'.
 // Note that Windows typically doesn't know how to open .ppm
 // images. Use Gimp or any other seious image processing
 // software to display .ppm images.
 // Also, note that because of Windows file format management,
 // you may have to modify this file to get image output on
 // Windows machines to work properly.
 //  
  
 FILE *f;
 unsigned char *bits24;
 double *rgbIm;
 
 if (im!=NULL)
  if (im->rgbdata!=NULL)
  {
   rgbIm=(double *)im->rgbdata;
   bits24=(unsigned char *)calloc(im->sx*im->sy*3,sizeof(unsigned char));
   for (int i=0; i<im->sx*im->sy*3; i++)
     *(bits24+i)=(unsigned char)(255.0*(*(rgbIm+i)));
   f=fopen(filename,"wb+");
   if (f==NULL)
   {
    fprintf(stderr,"Unable to open file %s for output! No image written\n",filename);
    return;
   }
   fprintf(f,"P6\n");
   fprintf(f,"# Output from PathTracer.c\n");
   fprintf(f,"%d %d\n",im->sx,im->sy);
   fprintf(f,"255\n");
   fwrite(bits24,im->sx*im->sy*3*sizeof(unsigned char),1,f);
   fclose(f);
   return;
   
   free(bits24);
  }
 fprintf(stderr,"imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{
 // De-allocates memory reserved for the image stored in 'im'
 if (im!=NULL)
 {
  if (im->rgbdata!=NULL) free(im->rgbdata);
  free(im);
 }
}

void dataOutput(double *im, int sx, char *name)
{
 FILE *f;
 double *imT;
 double HDRhist[1000];
 int i,j;
 double mx,mi,biw,pct;
 unsigned char *bits24;
 char pfmname[1024];
 
 imT=(double *)calloc(sx*sx*3,sizeof(double));
 memcpy(imT,im,sx*sx*3*sizeof(double));
 strcpy(&pfmname[0],name);
 strcat(&pfmname[0],".pfm");
 
  // Output the floating point data so we can post-process externally
 f=fopen(pfmname,"w");
 fprintf(f,"PF\n");
 fprintf(f,"%d %d\n",sx,sx);
 fprintf(f,"%1.1f\n",-1.0);
 fwrite(imT,sx*sx*3*sizeof(double),1,f);
 fclose(f);

 // Post processing HDR map - find reasonable cutoffs for normalization
 for (j=0; j<1000; j++) HDRhist[j]=0;

 mi=10e6;
 mx=-10e6;
 for (i=0; i<sx*sx*3; i++)
 {
  if (*(imT+i)<mi) mi=*(imT+i);
  if (*(imT+i)>mx) mx=*(imT+i);
 }

 for (i=0; i<sx*sx*3; i++)
 {
  *(imT+i)=*(imT+i)-mi;
  *(imT+i)=*(imT+i)/(mx-mi);
 }
 fprintf(stderr,"Image stats: Minimum=%f, maximum=%f\n",mi,mx);
 biw=1.000001/1000.0;
 // Histogram
 for (i=0; i<sx*sx*3; i++)
 {
  for (j=0;j<1000; j++)
   if (*(imT+i)>=(biw*j)&&*(imT+i)<(biw*(j+1))) {HDRhist[j]++; break;}
 }

 pct=.005*(sx*sx*3);
 mx=0;
 for (j=5; j<990;j++)
 {
  mx+=HDRhist[j];
  if (HDRhist[j+5]-HDRhist[j-5]>pct) break;
  if (mx>pct) break;
 }
 mi=(biw*(.90*j));

 for (j=990; j>5; j--)
 {
  if (HDRhist[j-5]-HDRhist[j+5]>pct) break;
 }
 mx=(biw*(j+(.25*(999-j))));

 fprintf(stderr,"Limit values chosen at min=%f, max=%f... normalizing image\n",mi,mx);

 for (i=0; i<sx*sx*3; i++)
 {
  *(imT+i)=*(imT+i)-mi;
  *(imT+i)=*(imT+i)/(mx-mi);
  if (*(imT+i)<0.0) *(imT+i)=0.0;
  if (*(imT+i)>1.0) *(imT+i)=1.0;
  *(imT+i)=pow(*(imT+i),.75);
 }

 bits24=(unsigned char *)calloc(sx*sx*3,sizeof(unsigned char));
 for (int i=0; i<sx*sx*3; i++)
  *(bits24+i)=(unsigned char)(255.0*(*(imT+i)));
 f=fopen(name,"wb+");
 if (f==NULL)
 {
  fprintf(stderr,"Unable to open file %s for output! No image written\n",name);
  return;
 }
 fprintf(f,"P6\n");
 fprintf(f,"# Output from PathTracer.c\n");
 fprintf(f,"%d %d\n",sx,sx);
 fprintf(f,"255\n");
 fwrite(bits24,sx*sx*3*sizeof(unsigned char),1,f);
 fclose(f);
 return;
   
 free(bits24);
 free(imT);
 
}

void cleanup(struct object3D *o_list, struct textureNode *t_list)
{
 // De-allocates memory reserved for the object list and for any loaded textures
 // Note that *YOU* must de-allocate any memory reserved for images
 // rendered by the raytracer.
 struct object3D *p, *q;
 struct textureNode *t, *u;

 p=o_list;		// De-allocate all memory from objects in the list
 while(p!=NULL)
 {
  q=p->next;
  free(p);
  p=q;
 }

 t=t_list;		// Delete texture Images
 while(t!=NULL)
 {
  u=t->next;
  if (t->im->rgbdata!=NULL) free(t->im->rgbdata);
  free(t->im);
  free(t);
  t=u;
 }
}

