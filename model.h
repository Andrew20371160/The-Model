
#ifndef _model_h_
#define _model_h_

#include <windows.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <cuda_runtime.h>
#define ROWS 720

#define COLS 1280
#define M_PI 3.14159
//threshold for differenece between angles
const FLOAT degree_threshold = 0.087f * 2;

//threshold for small forces
//used in insertion
const FLOAT force_threshold = 3;

//threshold for error in pixels count
const INT32 limit = 12;

const FLOAT zero_threshold = 0.001;
//threshold for small forces
//used in insertion

//forces vector size
const INT32 max_forces_size = 60;
const INT32 max_model_size = 30;
//threshold for small lengths
const FLOAT length_threshold = 10;

enum regions
{
    A = 2, B = 4, C = 6, D = 8
};

enum {
    DEFORM = 0, INSERT = 1, UNSTABLE = 2, STABLE = 3, REINIT
};
class point {

public:
    //coordinates
    INT32 x;
    INT32 y;
    //initializes the point to passed x,y
        // Default constructor
    __host__ __device__ point() : x(0), y(0) {}

    // Parameterized constructor
    __host__ __device__ point(const INT32 _x, const INT32 _y) : x(_x), y(_y) {}

    // Copy constructor
    __host__ __device__ point(const point& src_pt) : x(src_pt.x), y(src_pt.y) {}

    // Assignment operator
    __host__ __device__ point& operator=(const point& src_pt) {
        if (this != &src_pt) {
            x = src_pt.x;
            y = src_pt.y;
        }
        return *this;
    }
    __host__ __device__ point operator+(const point& src_pt) const {
        return point(x + src_pt.x, y + src_pt.y);
    }
    __host__ __device__ point operator-(const point& src_pt) const {
        return point(x - src_pt.x, y - src_pt.y);
    }
    __host__ __device__ INT32 dot(const point& src_pt)const {
        return x * src_pt.x + y * src_pt.y;
    }


    // Equality operator
    __host__ __device__ BOOL operator==(const point src_pt) const {
        return (x == src_pt.x) && (y == src_pt.y);
    }

    // Distance between two points
    __host__ __device__ FLOAT length(const point src_pt) const {
        return sqrtf((x - src_pt.x) * (x - src_pt.x) + (y - src_pt.y) * (y - src_pt.y));
    }

    // Midpoint between two points
    __host__ __device__ point mid_point(const point src_pt) const {
        return point((x + src_pt.x) / 2, (y + src_pt.y) / 2);
    }
    __host__ __device__ FLOAT normal_theta(const point src_pt)const;
    //returns theta in radians of the line connecting phantom & src_pt
    __host__ __device__ FLOAT theta(const point src_pt)const;
    //return dy/dx of the line from phantom to src_pt
    //midpoint between phantom and src_pt    // Distance between two points
    __host__ __device__ INT32 length_squared(const point src_pt) const {
        return ((x - src_pt.x) * (x - src_pt.x) + (y - src_pt.y) * (y - src_pt.y));
    }
};
BOOL line_eq(const point& p1, const point& p2, FLOAT& a, FLOAT& b);

//swaps 2 points

class model {

private:
    /*
    number of elements in the array
    */
    UINT32 elements_size;
    //one dimensional array of vertices
    point* vec;
    //an array that will be filled with forces each time a calculation is made
    FLOAT* force_vec;
    FLOAT* temp_vec;
    //stores current centroid value 
    point current_centroid;
public:
    /*
    constructors and desructors
    */
    //allocates 30 points for the model
    model(void);
    //destructor
    ~model(void);

    //returns number of elements in the model
    UINT32 get_size(void)const;

    /*
    insertion and removal
    */

    //inserts a point at an index
    BOOL insert(const point&, UINT32);
    //removes point at wanted index
    BOOL remove(UINT32 wanted_index);
    /*
        accessors
    */
    //constant accessor
    point& operator[](const UINT32)const;
    //variable accessor
    point& operator[](const UINT32);
    void sort(void);
    /*
    nagi's paper mathematics
    */
    BOOL forces(INT32* histogram);

    //performs deformation of vertecies
    BOOL deform(INT32* histogram);

    //performs insertion algorithm mentioned in the paper (lack of vertecies)
    void insert_vertecies_at_position(INT32* histogram);

    //performs deletion algorithm in the paper (too much vertecies)
    BOOL remove_extra_vertices(void);
    BOOL fix_intersection(void);
    BOOL   intersection_check(const point& p0, const point& p1, const point& p2, const point& p3, point& inter_p) const;

    //removes vertex between 2 lines that when difference between anlges of these lines is minimum
    BOOL remove_verticies_in_middle(void);
    BOOL remove_spikes(void);
    //removes vertex
    BOOL remove_small_lengths(void);

    point centroid(void)const;
    point centroid_average(void)const;
    FLOAT area(void)const;
    INT32  get_state(INT32* hist);

    //reinitializes the model where pt is the center
    void reinit(const point pt);
    BOOL set_centroid(void);
    //this function moves model vertices to accomodate the fact that the coordinates
    //changes when the camera head moves
    void move(const point& displacement); 
};

#endif