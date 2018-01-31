/*
CSC D18 - RayTracer code.

Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
Freely distributable for adacemic purposes only.

Uses Tom F. El-Maraghi's code for computing inverse
matrices. You will need to compile together with
svdDynamic.c

You need to understand the code provided in
this file, the corresponding header file, and the
utils.c and utils.h files. Do not worry about
svdDynamic.c, we need it only to compute
inverse matrices.

You only need to modify or add code in sections
clearly marked "TO DO" - remember to check what
functionality is actually needed for the corresponding
assignment!

Last updated: Aug. 2017   - F.J.E.
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


#include "utils.h"  // <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;
int sx;    // Size of the raytraced image
int photon_n=100000;
int photon_k=0;

void buildScene(void) {
#include "buildscene.c"   // <-- Import the scene definition!
}

void hierarchical(double M1[4][4],double depth,double ra, double rd, double rs, double rg,
    double alpha, double r_index,double shiny){
    if(depth>0){      
        for (int i=0;i<2;i++){
            double R = (double)rand() / (double)RAND_MAX;
            double G = (double)rand() / (double)RAND_MAX;
            double B = (double)rand() / (double)RAND_MAX;
            //struct object3D *o=newSphere(ra,rd,rs,rg,R,G,B,alpha,r_index,shiny);
            struct object3D *o=newCyl(ra,rd,rs,rg,R,G,B,alpha,r_index,shiny);
            memcpy(o->T,M1,16*sizeof(double));
            invert(&o->T[0][0],&o->Tinv[0][0]);
            insertObject(o,&object_list);
            
            double local[4][4]={{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
            ScaleMat(local,.75,.75,.75);
            RotateZMat(local, PI/3*pow(-1,i));	
            TranslateMat(local, 1.3, 1.3,0);
            matMult(M1,local);  
            hierarchical(local, depth-1, ra,rd,rs,rg,alpha,r_index,shiny);    
       }
    }
}

void
rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b,
        struct colourRGB *col) {
// This function implements the shading model as described in lecture. It takes
// - A pointer to the first object intersected by the ray (to get the colour properties)
// - The coordinates of the intersection point (in world coordinates)
// - The normal at the point
// - The ray (needed to determine the reflection direction to use for the global component, as well as for
//   the Phong specular component)
// - The current racursion depth
// - The (a,b) texture coordinates (meaningless unless texture is enabled)
//
// Returns:
// - The colour for this ray (using the col pointer)
//

    struct colourRGB tmp_col;  // Accumulator for colour components
    double R, G, B;      // Colour for the object in R G and B

// This will hold the colour as we process all the components of
// the Phong illumination model
    tmp_col.R = 0;
    tmp_col.G = 0;
    tmp_col.B = 0;

    if (obj->texImg == NULL)   // Not textured, use object colour
    {
        R = obj->col.R;
        G = obj->col.G;
        B = obj->col.B;
    } else {
// Get object colour from the texture given the texture coordinates (a,b), and the texturing function
// for the object. Note that we will use textures also for Photon Mapping.
        obj->textureMap(obj->texImg, a, b, &R, &G, &B);
    }

//////////////////////////////////////////////////////////////
// TO DO: Implement this function. Refer to the notes for
// details about the shading model.
//////////////////////////////////////////////////////////////
    double hasLS = 0;

    struct object3D *currentLS = object_list;
    struct point3D *s = (struct point3D *) calloc(1, sizeof(struct point3D));
    while (currentLS) {
        if (currentLS->isLightSource) {
            hasLS = 1;
            int unblock = 0;

            int rand_k = 5;
            for (int k = 0; k < rand_k; k++) {
                double lsx, lsy, lsz;
                (currentLS->randomPoint)(currentLS,    &lsx, &lsy, &lsz);
                // intersection point to light source, not normalize direction
                s = newPoint(lsx - p->px, lsy - p->py, lsz - p->pz);

                double lambda_pl;
                double a1, b1;
                struct object3D *temp_obj;
                struct point3D p1, n1;
                struct ray3D *r_pl = newRay(p, s);
                findFirstHit(r_pl, &lambda_pl, obj, &temp_obj, &p1, &n1, &a1, &b1);
                free(r_pl);

                if (lambda_pl < 0 || lambda_pl > 1) {
                    unblock++;
                }
            }

            if (unblock > 0) {
                double ratio = (double) unblock / (double) rand_k;
                // Phong modele
                double Isd_R = R * currentLS->col.R * ratio;
                double Isd_G = G * currentLS->col.G * ratio;
                double Isd_B = B * currentLS->col.B * ratio;

                // Diffuse term
                normalize(s);
                double dot_ns = dot(n, s);
                if (obj->alb.rd > 0) {
                    double max_res1 = obj->frontAndBack == 1 || obj->alpha < 1 ? fabs(dot_ns) : max(0, dot_ns);
                    tmp_col.R += obj->alb.rd * Isd_R * max_res1;
                    tmp_col.G += obj->alb.rd * Isd_G * max_res1;
                    tmp_col.B += obj->alb.rd * Isd_B * max_res1;
                }

                // Specular hightlights term
                if (obj->alb.rs > 0) {
                    struct point3D *m = newPoint(n->px * 2 * dot_ns - s->px, n->py * 2 * dot_ns - s->py,
                                                 n->pz * 2 * dot_ns - s->pz);

                    // intersection point to camera center, norm
                    struct point3D *c = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
                    normalize(m);
                    double dot_cm = dot(c, m);
                    double max_res2 = max(0, dot_cm);
                    tmp_col.R += obj->alb.rs * Isd_R * pow(max_res2, obj->shinyness);
                    tmp_col.G += obj->alb.rs * Isd_G * pow(max_res2, obj->shinyness);
                    tmp_col.B += obj->alb.rs * Isd_B * pow(max_res2, obj->shinyness);

                    free(c);
                    free(m);
                }
            }
        }

        currentLS = currentLS->next;
    }
    free(s);

    if (obj->alphaMap) alphaMap(obj->alphaMap, a, b, &obj->alpha);


    if (hasLS) {
        obj->alb.ra = 0.3;

        //ambient term
        tmp_col.R += obj->alb.ra * R;
        tmp_col.G += obj->alb.ra * G;
        tmp_col.B += obj->alb.ra * B;

//semi transparent object and refraction rays should be implemented
        if (obj->alpha < 1) {
            tmp_col.R *= obj->alpha;
            tmp_col.G *= obj->alpha;
            tmp_col.B *= obj->alpha;

            double r_idx1 = ray->goInside ? 1.0 : obj->r_index;
            double r_idx2 = ray->goInside ? obj->r_index : 1.0;

            struct point3D *op_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
            double cos_theta1 = dot(op_d, n);
            double sin_theta1 = sqrt(1 - pow(cos_theta1, 2));
            double sin_theta2 = (double) (r_idx1 / r_idx2) * sin_theta1;

            //if not total internal reflection, refraction angle < 90
            if (sin_theta2 < 1 && sin_theta2 > 0) {
                double temp1 = r_idx1 / r_idx2;
                double temp2 = -dot(n, &ray->d);
                double temp3 = temp1 * temp2 - pow(1 - temp1 * temp1 * (1 - temp2 * temp2), 0.5);
                struct point3D *d_refract = newPoint(temp1 * ray->d.px + temp3 * n->px,
                                                     temp1 * ray->d.py + temp3 * n->py,
                                                     temp1 * ray->d.pz + temp3 * n->pz);
                normalize(d_refract);

                struct colourRGB refracted_color;
                struct ray3D *refracted_ray = newRay(p, d_refract);
                refracted_ray->goInside = 1 - ray->goInside;
                rayTrace(refracted_ray, depth + 1, &refracted_color, obj);

                tmp_col.R += (1 - obj->alpha) * R * refracted_color.R;
                tmp_col.G += (1 - obj->alpha) * B * refracted_color.G;
                tmp_col.B += (1 - obj->alpha) * G * refracted_color.B;

                free(d_refract);
                free(refracted_ray);
            }

            free(op_d);
        }

// global reflection component
// trace rays from the surface point to possible incident directions
// if obj has specular reflection
        if (obj->alb.rg > 0) {
            struct colourRGB new_col;
            double dot_rdn = dot(&ray->d, n);
            struct point3D *ms = newPoint(n->px * (-2) * dot_rdn + ray->d.px, n->py * (-2) * dot_rdn + ray->d.py,
                                          n->pz * (-2) * dot_rdn + ray->d.pz);
            normalize(ms);
            struct ray3D *new_ray = newRay(p, ms);
            rayTrace(new_ray, depth + 1, &new_col, obj);

            tmp_col.R += obj->alb.rg * new_col.R;
            tmp_col.G += obj->alb.rg * new_col.G;
            tmp_col.B += obj->alb.rg * new_col.B;

            free(ms);
            free(new_ray);
        }
    }

    // apply the photon colour to current colour
    if (obj->photonMap) {
        double photon_R, photon_G, photon_B;
        obj->textureMap(obj->photonMap, a, b, &photon_R, &photon_G, &photon_B);

        tmp_col.R += obj->alb.rd*(obj->alpha)*photon_R*((double)photon_k/photon_n);
        tmp_col.G += obj->alb.rd*(obj->alpha)*photon_G*((double)photon_k/photon_n);
        tmp_col.B += obj->alb.rd*(obj->alpha)*photon_B*((double)photon_k/photon_n);
	}

// Be sure to update 'col' with the final colour computed here!
    col->R = tmp_col.R > 1 ? 1 : tmp_col.R;
    col->G = tmp_col.G > 1 ? 1 : tmp_col.G;
    col->B = tmp_col.B > 1 ? 1 : tmp_col.B;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p,
                  struct point3D *n, double *a, double *b) {
// Find the closest intersection between the ray and any objects in the scene.
// Inputs:
//   *ray    -  A pointer to the ray being traced
//   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
//              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
//              projection
// Outputs:
//   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
//   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
//              this ray (this is required so you can do the shading)
//   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
//   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
//   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

/////////////////////////////////////////////////////////////
// TO DO: Implement this function. See the notes for
// reference of what to do in here
/////////////////////////////////////////////////////////////
    double lambda_res = -1.0;
    *lambda = lambda_res;
    struct point3D p_res;
    struct point3D n_res;
    double a_res, b_res;

    struct object3D *currentObj = object_list;
    while (currentObj) {
// ignore self-intersections and for rays originating at the center of projection
        if ((currentObj != Os) && !currentObj->isLightSource) {
            (currentObj->intersect)(currentObj, ray, &lambda_res, &p_res, &n_res, &a_res, &b_res);

            // get the smaller lambda
            if (lambda_res > 0 && (lambda_res < *lambda || *lambda <= 0)) {
                *lambda = lambda_res;
                memcpy(p, &p_res, sizeof(struct point3D));
                memcpy(n, &n_res, sizeof(struct point3D));
                *a = a_res;
                *b = b_res;
                *obj = currentObj;
            }
        }
        currentObj = currentObj->next;
    }
}

void tracePhotonRay(int depth, struct ray3D *ray, struct object3D *os, double R, double G, double B) {
    if (depth > MAX_DEPTH) {
        return;
    }

    double lambda_r, ar, br;
    struct object3D *hit_obj;
    struct point3D pr, nr;
    findFirstHit(ray, &lambda_r, os, &hit_obj, &pr, &nr, &ar, &br);

    if (lambda_r > 0) {
        if (hit_obj->alpha < 1) {
            // if ray hits a refracting object
            double r_idx1 = ray->goInside ? 1.0 : hit_obj->r_index;
            double r_idx2 = ray->goInside ? hit_obj->r_index : 1.0;
            struct point3D *op_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
            double cos_theta1 = dot(op_d, &nr);
            double sin_theta1 = sqrt(1 - pow(cos_theta1, 2));
            double sin_theta2 = (double) (r_idx1 / r_idx2) * sin_theta1;

            if (sin_theta2 < 1 && sin_theta2 > 0) {
                double temp1 = r_idx1 / r_idx2;
                double temp2 = -dot(&nr, &ray->d);
                double temp3 = temp1 * temp2 - pow(1 - temp1 * temp1 * (1 - temp2 * temp2), 0.5);
                struct point3D *d_refract = newPoint(temp1 * ray->d.px + temp3 * nr.px,
                                                     temp1 * ray->d.py + temp3 * nr.py,
                                                     temp1 * ray->d.pz + temp3 * nr.pz);
                normalize(d_refract);
                struct ray3D *refract_ray = newRay(&pr, d_refract);
                refract_ray->goInside = 1 - ray->goInside;

                double alp = 1 - hit_obj->alpha;
                tracePhotonRay(depth+1, refract_ray, hit_obj, alp * hit_obj->col.R * R, alp * hit_obj->col.G * G,
                               alp * hit_obj->col.B * B);

                free(d_refract);
                free(refract_ray);
            }

            free(op_d);

        } else if (hit_obj->alb.rg > 0) {
            // if ray hits a reflecting object
            double dot_rdn = dot(&ray->d, &nr);
            struct point3D *ms = newPoint(nr.px * (-2) * dot_rdn + ray->d.px, nr.py * (-2) * dot_rdn + ray->d.py,
                                          nr.pz * (-2) * dot_rdn + ray->d.pz);
            normalize(ms);
            ray = newRay(&pr, ms);
            tracePhotonRay(depth+1, ray, hit_obj, hit_obj->alb.rg * R, hit_obj->alb.rg * G, hit_obj->alb.rg * B);
            free(ms);
        }

        if (hit_obj->alb.rd && depth > 0) {
            // ray hits a diffuse surface & been reflected/refracted
            // store one photon at the location of the intersection point and stop tracing
            double *photon_rgb = (double *) hit_obj->photonMap->rgbdata;

            int i = ar * sx;
            int j = br * sx;

            photon_rgb[3 * (i + sx * j)] += R;
            photon_rgb[3 * (i + sx * j) + 1] += G;
            photon_rgb[3 * (i + sx * j) + 2] += B;

            photon_k =photon_k+1;
        }
    }
}

void setPixelCol(double *rgbIm, int r, int dx) {
    for (int j = 0; j < sx; j++) {
        for (int i = 0; i < sx; i++) {
            int c = -r;

            double col_r = rgbIm[3 * (i + sx * j)];
            double col_g = rgbIm[3 * (i + sx * j) + 1];
            double col_b = rgbIm[3 * (i + sx * j) + 2];

            while (c < r) {
                int x = dx == 1 ? i+c : i;
                x = x < 0 ? x+sx : x;
                x = x >= sx ? x-sx:x;

                int y = dx == 0 ? j+c : j;
                y = y < 0 ? y+sx : y;
                y = y >= sx ? y-sx:y;

                col_r += rgbIm[3 * (x + sx * y)];
                col_g += rgbIm[3 * (x + sx * y) + 1];
                col_b += rgbIm[3 * (x + sx * y) + 2];
                c++;
            }

            col_r = (double)col_r/(2*r+1);
            col_g = (double)col_g/(2*r+1);
            col_b = (double)col_b/(2*r+1);

            rgbIm[3 * (i + sx * j)] = col_r;
            rgbIm[3 * (i + sx * j) + 1] = col_g;
            rgbIm[3 * (i + sx * j) + 2] = col_b;
        }
    }
}

void boxBlur(int r) {
    // r is the radius of the photon circle
    // reference http://amritamaz.net/blog/understanding-box-blur
    r = abs(2*r) > sx ? sx/2 : r;

    struct object3D *currentObj = object_list;
    while (currentObj) {
        if (currentObj->photonMap) {
            double *rgbIm = (double *)currentObj->photonMap->rgbdata;

            // for every neighboring pixel within radius r in x direction
            setPixelCol(rgbIm, r, 1);
            // for every neighboring pixel within radius r in y direction
            setPixelCol(rgbIm, r, 0);
        }

        currentObj = currentObj->next;
    }
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os) {
// Trace one ray through the scene.
//
// Parameters:
//   *ray   -  A pointer to the ray being traced
//   depth  -  Current recursion depth for recursive raytracing
//   *col   - Pointer to an RGB colour structure so you can return the object colour
//            at the intersection point of this ray with the closest scene object.
//   *Os    - 'Object source' is a pointer to the object from which the ray
//            originates so you can discard self-intersections due to numerical
//            errors. NULL for rays originating from the center of projection.
    double lambda;   // Lambda at intersection
    double a, b;    // Texture coordinates
    struct object3D *obj;  // Pointer to object at intersection
    struct point3D p, n;
    struct colourRGB I;  // Colour returned by shading function

    if (depth > MAX_DEPTH) // Max recursion depth reached. Return invalid colour.
    {
        col->R = -1;
        col->G = -1;
        col->B = -1;
        return;
    }

///////////////////////////////////////////////////////
// TO DO: Complete this function. Refer to the notes
// if you are unsure what to do here.
///////////////////////////////////////////////////////
    findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

    if (lambda > 0) {
        ray->goInside = Os == obj ? 0 : 1;
        rtShade(obj, &p, &n, ray, depth, a, b, &I);
        col->R = max(0, I.R);
        col->G = max(0, I.G);
        col->B = max(0, I.B);
    } else {
// background colour
        col->R = 0;
        col->G = 0;
        col->B = 0;
    }

}

int main(int argc, char *argv[]) {
// Main function for the raytracer. Parses input parameters,
// sets up the initial blank image, and calls the functions
// that set up the scene and do the raytracing.
    struct image *im;  // Will hold the raytraced image
    struct view *cam;  // Camera and view for this scene
    int antialiasing;  // Flag to determine whether antialiaing is enabled or disabled
    char output_name[1024];  // Name of the output file for the raytraced .ppm image
    struct point3D e;    // Camera view parameters 'e', 'g', and 'up'
    struct point3D g;
    struct point3D up;
    double du, dv;     // Increase along u and v directions for pixel coordinates
    struct ray3D *ray;   // Structure to keep the ray from e to a pixel
    struct colourRGB col;    // Return colour for raytraced pixels
    struct colourRGB background;   // Background colour
    int i, j;     // Counters for pixel coordinates
    unsigned char *rgbIm;

    if (argc < 5) {
        fprintf(stderr, "RayTracer: Can not parse input parameters\n");
        fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
        fprintf(stderr, "   size = Image size (both along x and y)\n");
        fprintf(stderr, "   rec_depth = Recursion depth\n");
        fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
        fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }
    sx = atoi(argv[1]);
    MAX_DEPTH = atoi(argv[2]);
    if (atoi(argv[3]) == 0) antialiasing = 0; else antialiasing = 1;
    strcpy(&output_name[0], argv[4]);

    fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
    fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
    if (!antialiasing) fprintf(stderr, "Antialising is off\n");
    else fprintf(stderr, "Antialising is on\n");
    fprintf(stderr, "Output file name: %s\n", output_name);

    object_list = NULL;
    light_list = NULL;
    texture_list = NULL;

// Allocate memory for the new image
    im = newImage(sx, sx);
    if (!im) {
        fprintf(stderr, "Unable to allocate memory for raytraced image\n");
        exit(0);
    } else rgbIm = (unsigned char *) im->rgbdata;

///////////////////////////////////////////////////
// TO DO: You will need to implement several of the
//        functions below. For Assignment 2, you can use
//        the simple scene already provided. But
//        for Assignment 3 you need to create your own
//        *interesting* scene.
///////////////////////////////////////////////////
    buildScene();    // Create a scene. This defines all the
    // objects in the world of the raytracer
//////////////////////////////////////////
// TO DO: For Assignment 2 you can use the setup
//        already provided here. For Assignment 3
//        you may want to move the camera
//        and change the view parameters
//        to suit your scene.
//////////////////////////////////////////

// Mind the homogeneous coordinate w of all vectors below. DO NOT
// forget to set it to 1, or you'll get junk out of the
// geometric transformations later on.

// Camera center is at (0,0,-1)
    e.px = 0;
    e.py = 0;
    e.pz = -3;
    e.pw = 1;

// To define the gaze vector, we choose a point 'pc' in the scene that
// the camera is looking at, and do the vector subtraction pc-e.
// Here we set up the camera to be looking at the origin.
    g.px = 0 - e.px;
    g.py = 0 - e.py;
    g.pz = 0 - e.pz;
    g.pw = 1;
// In this case, the camera is looking along the world Z axis, so
// vector w should end up being [0, 0, -1]

// Define the 'up' vector to be the Y axis
    up.px = 0;
    up.py = 1;
    up.pz = 0;
    up.pw = 1;

// Set up view with given the above vectors, a 4x4 window,
// and a focal length of -1 (why? where is the image plane?)
// Note that the top-left corner of the window is at (-2, 2)
// in camera coordinates.
    cam = setupView(&e, &g, &up, -1, -2, 2, 4);

    if (cam == NULL) {
        fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
        cleanup(object_list, light_list, texture_list);
        deleteImage(im);
        exit(0);
    }

// Set up background colour here
    background.R = 0;
    background.G = 0;
    background.B = 0;

// Do the raytracing
//////////////////////////////////////////////////////
// TO DO: You will need code here to do the raytracing
//        for each pixel in the image. Refer to the
//        lecture notes, in particular, to the
//        raytracing pseudocode, for details on what
//        to do here. Make sure you undersand the
//        overall procedure of raytracing for a single
//        pixel.
//////////////////////////////////////////////////////
    du = cam->wsize / (sx - 1);    // du and dv. In the notes in terms of wl and wr, wt and wb,
    dv = -cam->wsize / (sx - 1);   // here we use wl, wt, and wsize. du=dv since the image is
    // and dv is negative since y increases downward in pixel
    // coordinates and upward in camera coordinates.

    fprintf(stderr, "View parameters:\n");
    fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
    fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
    printmatrix(cam->C2W);
    fprintf(stderr, "World to camera conversion matrix:\n");
    printmatrix(cam->W2C);
    fprintf(stderr, "\n");



    //photon mapping
    fprintf(stderr, "Creating photon map... Number of Random Ray: %d \n", photon_n);

    // initial photon map for each object
    struct object3D *currentObj = object_list;
    while (currentObj) {
        if (!currentObj->isLightSource && (currentObj->alb.rd > 0)) {
            currentObj->photonMap = (struct image *) calloc(1, sizeof(struct image));
            if (currentObj->photonMap != NULL) {
                currentObj->photonMap->rgbdata = NULL;
                currentObj->photonMap->sx = sx;
                currentObj->photonMap->sy = sx;
                currentObj->photonMap->rgbdata = (double *) calloc(sx * sx * 3, sizeof(double));
            }
        }
        currentObj = currentObj->next;
    }

    // find photon center
    currentObj = object_list;
    int ls_n = 0;
    while (currentObj) {
        if (currentObj->isLightSource) {
            ls_n ++;
            for (int rk = 0; rk < photon_n; rk++) {
                // random light with random direction
                struct ray3D rand_ray;
                (currentObj->lsRandomRay)(currentObj, &rand_ray, &e);

               tracePhotonRay(1, &rand_ray, currentObj, currentObj->col.R, currentObj->col.G, currentObj->col.B);
            }
        }
        currentObj = currentObj->next;
    }
    photon_n *= ls_n;

    // generate photons with a sphere
    // loop several times to make the photon more smooths
    for(int bn = 0; bn < 3; bn ++)
    {
        boxBlur(7);
    }
  

    // For each of the pixels in the image
    fprintf(stderr, "\nRendering row: ");

    struct point3D *pc;
    struct point3D *d;
    struct colourRGB new_col;
    int anti_k = 5;
    
#pragma omp parallel for schedule(dynamic, 32) shared(anti_k, rgbIm, object_list) private(j)
    for (j = 0; j < sx; j++)   // For each of the pixels in the image
    {
         fprintf(stderr, "%d, ", j);
        //fprintf(stderr, ".");
#pragma omp parallel for private(pc, d, ray, col, new_col, i)

        for (i = 0; i < sx; i++) {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing!
            ///////////////////////////////////////////////////////////////////
            // anti-aliasing, random number within the pixel square

            col.R = 0;
            col.B = 0;
            col.G = 0;

            for (int ak = 0; ak < anti_k; ak++) {
                // get random x and y within range
                double pc_x = du * (double) rand() / RAND_MAX + ((-cam->wsize / 2) + i * du - du / 2);
                double pc_y = dv * (double) rand() / RAND_MAX + ((cam->wsize / 2) + j * dv - dv / 2);
                pc = newPoint(pc_x, pc_y, cam->f);
                matVecMult(cam->C2W, pc); // pij in world coordinate

                d = newPoint(pc->px - cam->e.px, pc->py - cam->e.py, pc->pz - cam->e.pz);
                normalize(d);

                ray = newRay(&cam->e, d);


                rayTrace(ray, 1, &new_col, NULL);
                col.R += new_col.R;
                col.G += new_col.G;
                col.B += new_col.B;

                free(ray);
                free(pc);
                free(d);
            }

            col.R = col.R / anti_k;
            col.G = col.G / anti_k;
            col.B = col.B / anti_k;

            rgbIm[3 * (i + sx * j)] = (unsigned char) (col.R * 255);
            rgbIm[3 * (i + sx * j) + 1] = (unsigned char) (col.G * 255);
            rgbIm[3 * (i + sx * j) + 2] = (unsigned char) (col.B * 255);
        } // end for i
    } // end for j

    fprintf(stderr, "\nDone!\n");

// Output rendered image
    imageOutput(im, output_name);

// Exit section. Clean up and return.
    cleanup(object_list, light_list, texture_list);    // Object, light, and texture lists
    deleteImage(im);         // Rendered image
    free(cam);           // camera view
    exit(0);
}
