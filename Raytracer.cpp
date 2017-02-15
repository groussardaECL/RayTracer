// Cours RayTracer http://goo.gl/40Etg5
// Utilisation de const https://openclassrooms.com/courses/la-const-correctness-expliquee-aux-zeros
// Utilisation du pointer this http://en.cppreference.com/w/cpp/language/this
// Pour scene et le nombre de spheres variables http://stackoverflow.com/questions/1657883/variable-number-of-arguments-in-c
// Utilisation de Data Structure http://www.cplusplus.com/doc/tutorial/structures/
// Obj Charger http://pastebin.com/DHftXgb5

// Includes libraries 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "CImg.h"
#include <stdarg.h>
#include <random>
#include<ctime>
std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0, 1);



using namespace std;
#define PI 3.14159265359

////////////////////////////////////////////////////////////
//Class Declarations////////////////////////////////////////
////////////////////////////////////////////////////////////

class Vector {
public:
	Vector(double x, double y, double z);
	Vector(const Vector &b);
	const double operator[](int i) const;
	double squaredNorm() const;
	void Normalize();
private:
	double xyz[3];
};

class Camera {
public:
	Camera(const Vector &origin_cam, int width_cam, int height_cam, double aov_d_cam);
	Vector origin;
	int W;
	int H;
	double aov_d;
	double aov;
	//vector<unsigned char> pixels(int n,int init); See definition
};

class Ray {
public:
	Ray(const Vector &origin_ray, const Vector &direction_ray);
	Ray(const Camera &camera, int x, int y); // Special builder for the rays out of the camera (only in (0,0,-1) direction)
	Vector origin;
	Vector direction; //Normalized Vector of direction_ray
};

class Material {
public:
	Material();
	Material(double R_mat, double G_mat, double B_mat);
	Material(int property, double R_mat, double G_mat, double B_mat);
	Material(double coeff_d, double coeff_r, double coeff_t, double R_mat, double G_mat, double B_mat);
	Vector color;
	double coeff_diffuse;
	double coeff_reflective;
	double coeff_transparent;
	double coeff_refraction;
};

class Sphere {
public:
	Sphere(const Vector &origin_sphr, int radius_sphr, const Material &M_sphr);
	double getIntersection(const Ray &R); // Returns data structure resultIntersection with index, P, T, and bool.
	Vector origin;
	int radius;
	Material M;
};

class Light {
public:
	Light(const Vector &origin_light, double intensity_light);
	Vector origin;
	double intensity;
};

struct dataIntersect {
	int index;
	double t;
	bool hasIntersection;
	//Vector P; Pour plus tard
};

class Scene {
public:
	Scene(vector<Sphere>,vector<Light>);
	vector<Sphere> VectSphr;
	vector<Light> VectLght;
	int CountSphr();
	int CountLght();
	dataIntersect getSceneIntersection(const Ray &R);
	Vector getColor(const Ray &R, int max_bounces);
	double coeff_refraction_scene = 1.;
};

// Operators Overload
Vector operator+(const Vector &a, const Vector &b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
};
Vector operator-(const Vector &a, const Vector &b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
};
Vector operator*(const Vector &a, const Vector &b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
};
// Scalar product
double dot(const Vector &a, const Vector &b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};
// Product with a double
Vector operator*(double a, const Vector &b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector &b, double a) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
// Division with a double
Vector operator/(const Vector &b, double a) {
	return Vector(b[0] / a, b[1] / a, b[2] / a);
}
// Cross Product (wedge)
Vector cross(const Vector &a, const Vector &b) {
	return Vector(a[1] * b[2] - b[1] * a[2], a[2] * b[0] - b[2] * a[0], a[0] * b[1] - b[0] * a[1]);
}


////////////////////////////////////////////////////////////
//Main Program//////////////////////////////////////////////
////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) {
	time_t tstart;
	tstart = time(0);
	Camera cam(Vector(0, 0, 55), 1024, 1024, 80);
	vector<unsigned char> pixels(cam.W*cam.H * 3, 0);
	Sphere SphrBckGrnd	(Vector(0, 0, -1000),	940,	Material(0, 0., 255., 0.));
	Sphere SphrGrnd		(Vector(0, -1000, 0),	990,	Material(0, 0., 0., 255.));
	Sphere SphrTop		(Vector(0, 1000, 0),	940,	Material(0, 255., 255., 0.));
	Sphere SphrBck		(Vector(0, 0, 1000),	940,	Material(0, 255., 0., 255.));
	Sphere SphrRght		(Vector(1000, 0, 0),	940,	Material(0, 255., 0., 0.));
	Sphere SphrLft		(Vector(-1000, 0, 0),	940,	Material(0, 0., 255., 255.));
	Sphere Sphr1		(Vector(0, 0, 10),		10,		Material(2, 255., 255., 255.));
	// Sphr1			(Vector(0, 0, 0),		10,		Material(50, 50, 0, 255., 255., 255.));
	//Sphere Sphr2		(Vector(0, 5, -30),		5,		Material(0, 255., 255., 0.));
	//Sphere Sphr3		(Vector(-10, 0, 40),	5,		Material(0, 255., 0., 255.));
	/*Sphere Sphr4		(Vector(-30, 30, 10),	5,		Material(1, 255., 255., 255.));
	Sphere Sphr5		(Vector(0, 30, 10),		5,		Material(1, 255., 255., 255.));
	Sphere Sphr6		(Vector(30, 30, 10),	5,		Material(1, 255., 255., 255.));
	Sphere Sphr7		(Vector(-30, -30, 10),	5,		Material(1, 255., 255., 255));
	Sphere Sphr8		(Vector(0, -30, 10),	5,		Material(1, 255., 255., 255.));
	Sphere Sphr9		(Vector(30, -30, 10),	5,		Material(1, 255., 255., 255.));*/
	Light Lght1			(Vector(-10, 20, 40),	255.*1200.);
	//Light Lght2			(Vector(10, 20, -10),	255.*300.);
	Scene Scn({ SphrBckGrnd,SphrGrnd,SphrTop,SphrBck,SphrLft,SphrRght,Sphr1 }, { Lght1 });
	bool display_true_intensity = false;
	double gamma = 2.2;
	int nbRaysPerPixel = 1;
	int nbMaxRebounds = 10;
	
#pragma omp parallel for // For parallel computation
	double sumETA = 0.;
	double ETA = 0.;
	for (int i = 0; i<cam.H; i++) {
		
		double perc = (double)i / (double)cam.H;
		double timer = difftime(time (0), tstart);
		if (i > 10)
		{
			sumETA += (timer / perc);
			ETA = sumETA / (i-10);
		}
		
		printf("Image : %i/%i (%i%%) - Temps/ETA : %i/%i\n", i, cam.H, (int) (perc * 100), (int) timer, (int) ETA);
		for (int j = 0; j<cam.W; j++) {
			
			if (i == 644 && j == 432) {
			printf("STOP");
			}
			Vector colorPix(0.,0.,0.);
			for (int l = 0; l < nbRaysPerPixel; l++) {
				Vector colorPix_tmp = Scn.getColor(Ray(cam, i, j), nbMaxRebounds);
				colorPix = colorPix + colorPix_tmp;
			}
			colorPix = colorPix / nbRaysPerPixel;

			double R = 0;
			double G = 0;
			double B = 0;

			if (!display_true_intensity)
			{
				R = min(255., colorPix[0]);
				G = min(255., colorPix[1]);
				B = min(255., colorPix[2]);
			}
			else
			{
				R = colorPix[0];
				G = colorPix[1];
				B = colorPix[2];
			}

			pixels[i*cam.W + j] = (unsigned char) (255. * pow(R/255.,1/gamma));
			pixels[i*cam.W + j + cam.W*cam.H] = (unsigned char) (255. * pow(G / 255., 1 / gamma));
			pixels[i*cam.W + j + 2 * cam.W*cam.H] = (unsigned char) (255. * pow(B / 255., 1 / gamma));
		}
	}
	cimg_library::CImg<unsigned char> cimg(&pixels[0], cam.W, cam.H, 1, 3);
	cimg.save("fichier57.bmp");
}

////////////////////////////////////////////////////////////
//Class Definitions/////////////////////////////////////////
////////////////////////////////////////////////////////////

Vector::Vector(double x = 0. , double y = 0., double z = 0.) { //Builder with default values
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
};

Vector::Vector(const Vector &b) {  //Another builder from an existing vector - Reference (const myClass &myObject) to the vector to avoid duplication in memory
	memcpy(xyz, b.xyz, 3 * sizeof(double)); //Copy in memory, with memory allocation
};

const double Vector::operator[](int i) const { //Method that return the cordinates with brackets
	return xyz[i];
};

/* Squared Euclidean Norm */
double Vector::squaredNorm() const {
	return xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
};

void Vector::Normalize() {
	double n = sqrt((*this).squaredNorm()); // this is a explicit pointer to the Vector being normalized - Useful for disambiguation
	xyz[0] /= n;
	xyz[1] /= n;
	xyz[2] /= n;
}

Camera::Camera(const Vector &origin_cam, int width_cam, int height_cam, double aov_d_cam) {
	origin = origin_cam;
	W = width_cam;
	H = height_cam;
	aov_d = aov_d_cam;
	// Conversion into radians, "fov" = Field of view in rad
	aov = aov_d*3.14 / 180;
	// vector<unsigned char> pixels(W*H * 3, 0); //Test in order to include vector in Camera
	//Errors : 
	//1>c:\users\groussard\documents\visual studio 2015\projects\raytracer\raytracer\raytracer.cpp(152): error C3867: 'Camera::pixels': non-standard syntax; use '&' to create a pointer to member
	//1>c:\users\groussard\documents\visual studio 2015\projects\raytracer\raytracer\raytracer.cpp(152) : error C2109 : subscript requires array or pointer type
};

Material::Material(void) {
	this->color = Vector(1.,1.,1.);
	this->coeff_diffuse = 1;
	this->coeff_reflective = 0;
	this->coeff_transparent = 0;
	this->coeff_refraction = 1.52;

}

Material::Material(double R_mat, double G_mat, double B_mat) {
	this->color = Vector(R_mat / 255., G_mat / 255., B_mat / 255.);
	this->coeff_diffuse = 1;
	this->coeff_reflective = 0;
	this->coeff_transparent = 0;
	this->coeff_refraction = 1.52;
};

Material::Material(int property, double R_mat, double G_mat, double B_mat) {
	this->color = Vector(R_mat / 255, G_mat / 255, B_mat / 255);
	this->coeff_refraction = 1.52;
	switch (property)
	{
	case 0 :
		coeff_diffuse = 1;
		coeff_reflective = 0;
		coeff_transparent = 0;
		break;

	case 1:
		coeff_diffuse = 0;
		coeff_reflective = 1;
		coeff_transparent = 0;
		break;

	case 2:
		coeff_diffuse = 0;
		coeff_reflective = 0;
		coeff_transparent = 1;
		break;

	default :
		coeff_diffuse = 1;
		coeff_reflective = 0;
		coeff_transparent = 0;
		break;
	}
};

Material::Material(double coeff_d, double coeff_r, double coeff_t, double R_mat, double G_mat, double B_mat) {
	
	this->color = Vector(R_mat / 255., G_mat / 255., B_mat / 255.);
	double sum = coeff_r + coeff_d + coeff_t + 1;
	this->coeff_diffuse = coeff_d / sum;
	this->coeff_reflective = coeff_r / sum;
	this->coeff_transparent = coeff_t / sum;
	this->coeff_refraction = 1.52;
};



Sphere::Sphere(const Vector &origin_sphr, int radius_sphr, const Material &M_sphr) {
	this->origin = origin_sphr;
	this->radius = radius_sphr;
	this->M = M_sphr;
}

double Sphere::getIntersection(const Ray &R) {
	// The intersection method of Sphere return t, the solution of the quadratic equation. 
	// If t = -1., there is no intersection with the sphere. If t<=0, the intersection is P = R.origin + t*R.direction.
	double a = 1.;
	double b = 2.*dot(R.direction, R.origin - this->origin);
	double c = (R.origin - this->origin).squaredNorm() - pow(this->radius, 2);
	double delta = b * b - 4 * a * c;
	if (delta >= 0.) { 
		// If delta < 0, no intersection
		// The second solution is mathematically greater than the first one. 
		//With the following boolean, we choose to use only this one, or to compare and then choose the shortest (in case that we didn't think of)
		bool uniqueSolution = false; // Change to false to use the 2 solutions and compare the one closest to the origin of the ray

		if (uniqueSolution)
		{
			double t = ((-b - sqrt(delta)) / (2 * a));
			return t;
		}
		else
		{
			double t1 = (-b - sqrt(delta)) / (2 * a); // First solution
			double t2 = (-b + sqrt(delta)) / (2 * a); //Second solution

			if (t1 > 0 && t2 > 0)
			{
				if (t1 <= t2)
				{
					return t1;
				}
				else
				{
					return t2;
				}
			}
			else if (t1 > 0) { return t1; }
			else if (t2 > 0) { return t2; }
			
		}
	}
	else
	{
		double t = -1.;
		return t;
	}
}

Ray::Ray(const Vector &origin_ray, const Vector &direction_ray) {
	this->origin = origin_ray;
	this->direction = direction_ray;
	this->direction.Normalize();
}

Ray::Ray(const Camera &camera, int x, int y) {
	this->origin = camera.origin;
	this->direction = Vector(y + 0.5 - camera.W / 2., camera.H / 2. - x + 0.5, -camera.W / (2.*tan(camera.aov / 2.)));
	this->direction.Normalize();
}

Light::Light(const Vector &origin_light, double intensity_light) {
	origin = origin_light;
	intensity = intensity_light;
}

Scene::Scene(vector<Sphere> V_Sphere,vector<Light> V_Light) {
	this->VectSphr = V_Sphere;
	this->VectLght = V_Light;
}

int Scene::CountSphr(void) {
	return (int)this->VectSphr.size();
}

int Scene::CountLght(void) {
	return (int)this->VectLght.size();
}

 dataIntersect Scene::getSceneIntersection(const Ray &R) {
	dataIntersect result;
	int nSphr = this->CountSphr();
	result.hasIntersection = false;
	result.index = -1;
	result.t = -1.;
	//vector<double> T;
	for (int k = 0; k < nSphr; k++) {
		double tmp = this->VectSphr[k].getIntersection(R);
		if (tmp > 0) {
			if (1 / tmp > 1 / result.t) {
				result.t = tmp;
				result.hasIntersection = true;
				result.index = k;
			}
		}
	}	
	return result;
}

 

 Vector Scene::getColor (const Ray &R, int max_bounces) {
	dataIntersect intersectData = this->getSceneIntersection(R); // Search of the intersection of the ray with elements of the scene
	// Initialization of the returned vector
	//printf("%i - ", max_bounces);
	Vector pixel(0., 0., 0.);
	
	if (intersectData.hasIntersection)
	{
		// Definition of the intersection point on the k-th sphere
		Vector P = R.origin + intersectData.t*R.direction;

		// Definition of the material of the point
		Material sphere_material = this->VectSphr[intersectData.index].M;
		
		// Definiton of Normal Unit Vector in P
		Vector n = (P - this->VectSphr[intersectData.index].origin);
		n.Normalize();

		// For reflective material
		if (sphere_material.coeff_reflective > 0. && max_bounces > 0)
		{
			// Definition of the Reflected Vector from P
			Vector r = R.direction - 2 * dot(R.direction, n) * n;
			Vector pixel_reflective = this->getColor(Ray(P + (1./1000.) * n, r), max_bounces-1);
			pixel = pixel + pixel_reflective * sphere_material.color * sphere_material.coeff_reflective;
		}

		if (sphere_material.coeff_transparent > 0 && max_bounces > 0)
		{
			double orientation = (dot(R.direction,n) < 0) - (dot(R.direction,n) > 0);

			Vector n_trans = orientation * n;
			double coeff_ratio = pow((this->coeff_refraction_scene / sphere_material.coeff_refraction), orientation);

			double var_test = 1 - pow(coeff_ratio,2) * (1 - pow(dot(R.direction, n_trans), 2));

			if (var_test > 0)
			{
				Vector ray_refrated_norm = coeff_ratio * (R.direction - dot(R.direction, n_trans) * n_trans);
				Vector ray_refracted_trans = -sqrt(1 - pow(coeff_ratio, 2) * (1 - pow(dot(R.direction, n_trans), 2))) * n_trans;
				Ray ray_exit(P - (1./1000.) * n_trans, ray_refracted_trans + ray_refrated_norm);

				Vector pixel_returned = this->getColor(ray_exit, max_bounces - 1);
				pixel = pixel + pixel_returned * sphere_material.color * sphere_material.coeff_transparent;
			}
			else
			{
				Vector r = R.direction - 2 * dot(R.direction, n_trans) * n_trans;
				
				Vector pixel_transparent = this->getColor(Ray(P + (1./1000.) * n, r), max_bounces - 1);
				pixel = pixel + pixel_transparent * sphere_material.color * sphere_material.coeff_transparent;
			}
		}
		
		// For diffuse material
		if (sphere_material.coeff_diffuse > 0.)
		{
			for (int k = 0; k < this->CountLght(); k++)
			{
				bool pixel_not_null = false;
				// Vector to Light
				Vector l = this->VectLght[k].origin - P;
				double d = sqrt(l.squaredNorm());
				l.Normalize();

				Ray rayP(P + (1./1000.) * n, l); // Definition of a ray from P to light L, to detect shadows
				dataIntersect intersectData_2 = this->getSceneIntersection(rayP); // Returns the closest intersection point to P
				
				if (intersectData_2.hasIntersection)
				{
					if (sqrt((intersectData_2.t * l).squaredNorm()) >= d) 
					{pixel_not_null = true;}
		
				}
				else  
				{pixel_not_null = true;}

				if (pixel_not_null)
				{
					double intensity = (this->VectLght[k].intensity * std::max(0., dot(n, l) / pow(d, 2)));
					pixel = pixel + sphere_material.color * intensity * sphere_material.coeff_diffuse;
				}
			}
			
			
			if (max_bounces > 0 )
			{
				// For indirect lights
				// Creation of a random ray
				double r1 = distrib(engine);
				double r2 = distrib(engine);
				Vector randVector(sqrt(1 - r2) * cos(2 * PI * r1), sqrt(1 - r2) * sin(2 * PI * r1), sqrt(r2)); // Already normalized, by definition.
				
				/*Vector tangent1;
				if (n[0] != 0 || n[1] != 0)
					tangent1 = Vector(n[1], n[0], 0);
				else
					tangent1 = Vector(0, n[2], -n[1]);
				tangent1.Normalize();*/


				// Construction of an orthonormed system v1, v2, v3
				Vector constructionVector(distrib(engine), distrib(engine), distrib(engine));
				constructionVector.Normalize();
				Vector n2(n); // Vector normal (~z) of the new base
				Vector n0 = cross(n2, constructionVector); // Vector (~x) of the new base
				Vector n1 = cross(n2, n0); // Vector (~y) of the new base
				n0.Normalize();
				n1.Normalize();

				// Expression of the randRay in the classic orthonormed system x, y, z
				Vector randRay_newBase(n0[0] * randVector[0] + n1[0] * randVector[1] + n2[0] * randVector[2], n0[1] * randVector[0] + n1[1] * randVector[1] + n2[1] * randVector[2], n0[2] * randVector[0] + n1[2] * randVector[1] + n2[2] * randVector[2]);
				Ray randRay(P + 1. / 1000. * n, randRay_newBase);
				Vector indirectColor = this->getColor(randRay, max_bounces - 1);
				pixel = pixel + indirectColor * max(0., dot(randRay.direction, n)) * sphere_material.color * sphere_material.coeff_diffuse / PI;
			}
		}
	}
	//if (max_bounces <= 0) { printf("\n END\n", max_bounces); }
	return pixel;
 }

