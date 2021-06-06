#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

//implementation of 3d Vector
struct vec {
	
	double x,y,z;
	vec(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	vec() {}
	vec(const vec &v2) {x = v2.x; y = v2.y; z = v2.z;}
	vec(vec* v2) { x=v2->x; y=v2->y; z=v2->z; }

	vec operator + (const vec o) {
		vec r;
		r.x = this->x + o.x;
		r.y = this->y + o.y;
		r.z = this->z + o.z;
		return r;
	}
	vec operator - (const vec o) {
		vec r;
		r.x = this->x - o.x;
		r.y = this->y - o.y;
		r.z = this->z - o.z;
		return r;
	}
	vec operator * (const vec o) {
		vec r;
		r.x = this->x * o.x;
		r.y = this->y * o.y;
		r.z = this->z * o.z;
		return r;
	}
	double mag(){
		return std::sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	}
	vec operator * (const double s){
		vec n;
		n.x = this->x*s;
		n.y = this->y*s;
		n.z = this->z*s;
		return n;
	}
	vec unit(){
		vec r;
		double len = mag();
		r.x = this->x/len;
		r.y = this->y/len;
		r.z = this->z/len;
		return r;
	}
	void step(vec direction, double t){
		vec dir = direction.unit();
		x += dir.x*t;
		y += dir.y*t;
		z += dir.z*t;
	}
	
	vec space_distortion(){ //This function will act on ray, tricking the object's max_dist and normal function into thinking the ray is in a different location
		vec space(this); //TODO: test if copy constructor is necessary when passing arguments to methods
		/* uncomment to repeat space
		while( space.z > 10 ){
			space.z-=10;
		}*/
		return space;
	}


};

//because I removed the color class due to the annoying integer math in favor of vec, this function returns a usable clipped color
int get_rgb(double c) {
	return (int) (std::max(0.0, std::min(255.0, c)));
}

//intended to include point lights and sun lights
struct light {
	vec is, id;
	vec loc;
	int type; //0 is point falloff, 1 is sun falloff
	light(vec _is, vec _id, vec _loc, int _type) : is(_is), id(_id), loc(_loc), type(_type){}
};

//phong shading variable storage
struct material {
	vec  ks, kd, ka;
	double a;
	material(vec specular_constant, vec diffuse_constant, vec ambient_constant, double specular_noise) :  ks(specular_constant), kd(diffuse_constant), ka(ambient_constant), a(specular_noise) {}
};


//generalized object containing material, location and TODO rotation information
class scene_object {
	public:
		//phong reflection model
		vec ks, kd, ka;
		double a;
		vec loc;
		scene_object() {}
		scene_object(vec _loc, material m) : loc(_loc) , kd(m.kd), ks(m.ks), ka(m.ka), a(m.a) {}
		scene_object(vec _loc,  vec _ks, vec _kd, vec _ka, double _a) : loc(_loc), ks(_ks), kd(_kd), ka(_ka), a(_a) {} 

		virtual double max_dist(vec ray){}
		virtual vec normal(vec ray){}
	};

//implementation of scene_object, with radius parameter and sphere distance function
class sphere: public scene_object{
	public:
		double rad;
		sphere(vec _loc, double _rad, material m) : scene_object(_loc, m), rad(_rad) {}
		sphere(vec _loc, double _rad, vec _ks, vec _kd, vec _ka, double _a) : scene_object(_loc,  _ks, _kd, _ka, _a), rad(_rad) {}
		double max_dist(vec ray) override {
			vec space = ray.space_distortion();
			return (loc-space).mag()-rad;
		}
		vec normal(vec ray) override {
			vec space = ray.space_distortion();
			return (space-loc).unit();
		}	

};

//simply a ground plane
class ground_plane: public scene_object{
	public:
		ground_plane(vec _loc, material m) : scene_object(_loc, m) {}
		ground_plane(vec _loc, vec _ks, vec _kd, vec _ka, double _a) : scene_object(_loc,  _ks, _kd, _ka, _a) {}
		double max_dist(vec ray) override {
			return ray.y-loc.y;
		}
		vec normal(vec ray) override {
			return vec(0, 1, 0);
		}
};

class mandlebulb: public scene_object{
	public:
		float Power;
		int Iterations;
		float Bailout;
		vec coloring;
		mandlebulb(vec _loc, material m, float _Power, vec _coloring) : scene_object(_loc, m),Power(_Power),Iterations(100),Bailout(100),coloring(_coloring) {}
		//this distance function is not my own. Credit: http://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/
		double max_dist(vec ray) override {;
		    ray = loc-ray;
			vec h = ray;
			float dr = 1.0;
			float r = 0.0;
			float min_orbit = 255;
			for (int i=0; i < Iterations; i++){
				r = h.mag();
				if(r>Bailout) break;
				float theta = acos(h.z/r);
				float phi = atan2(h.y,h.x);
				dr = pow(r,Power-1.0)*Power*dr+1.0;
				float zr = pow(r,Power);
				theta = theta*Power;
				phi=phi*Power;
				h = vec(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta))*zr;
				if(h.mag()<min_orbit) min_orbit = h.mag();
				h = h+ray;
			}
            //orbit trapping: not my idea, but my code
			ka = coloring*min_orbit;
			return 0.5*log(r)*r/dr;
		}
		vec normal(vec ray) override {
			//return (ray-loc).unit();
			return vec(0,1,0);
		}
};


//TODO Box shadows are currently broken
class box: public scene_object{
	public:
		double w,l,h;
		
		box(vec _loc, material m, vec dimm) : scene_object(_loc, m), w(dimm.x), h(dimm.y), l(dimm.z) {}

		double max_dist(vec ray) override {
			double x = std::abs(ray.x-loc.x);
			double y = std::abs(ray.y-loc.y);
			double z = std::abs(loc.z-ray.z); //left-handed coordinate system
			if(x<w && y<h && z<l){ return -1.0; } //check if inside first
			if(x<w){
				if(y<h){ return z-l; }
				if(z>l){ return std::sqrt(std::pow(y-h, 2) + std::pow(z-l,2)); }
				return y-h;
			}
			if(y<h){
				if(z>l){ return std::sqrt(std::pow(x-w,2)+std::pow(z-l,2)); }
				return x-w;
			}
			if(z<l){
				return std::sqrt(std::pow(x-w,2) + std::pow(y-h,2));
			}
			return (vec(x,y,z)-vec(w,h,l)).mag();
		}
		vec normal(vec ray) override {
			double x = std::abs(ray.x-loc.x);
			double y = std::abs(ray.y-loc.y);
			double z = std::abs(loc.z-ray.z); //left-handed coordinate system
			vec res(1,1,1);
			if(ray.x<loc.x) { res = res + vec(-2,0,0); }
			if(ray.y<loc.y) { res = res + vec(0,-2,0); }
			if(ray.z>loc.z) { res = res + vec(0,0,-2); }
			if(x<w){
				if(y<h){ return vec(0,0,-1)*res; }
				if(z>l){ return vec(0,h,-l).unit()*res; }
				return vec(0,1,0)*res;
			}
			if(y<h){
				if(z>l){ return vec(w,0,-l).unit()*res; }
				return vec(1,0,0)*res;
			}
			if(z<l){
				return vec(w,h,0).unit()*res;
			}
			return (vec(x,y,z)-vec(w,h,l)).unit()*res;

		}
};

//the dot product of two vectors
double dot(vec a, vec b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

vec WHITE(255,255,255);
vec BLACK(0,0,0);
vec GREEN(0,255,0);
vec BLUE(0,0,255);
vec RED(255,0,0);
vec GRAY = WHITE*.5;

//points to the right
vec I = vec(1,0,0);
//points up
vec J = vec(0,1,0);
//points into the screen plane
vec K = vec(0,0,1);

const double MIN_DIST = .001;
const double MAX_DIST = 800;
const int MAX_NUM_STEPS = 1000;
const double MAX_STEP = 100;
const int NUM_OBJECTS = 2;
const int NUM_LIGHTS = 2;
const double SHADOW_CORRECTION = .05;

const int W = 640;
const int H = 480;

const vec ORIGIN(0,0,0);

const double FOCAL_LENGTH = 1.0;
const double FRUSTUM_HEIGHT = H/300.0;
const double FRUSTUM_WIDTH = W/300.0;

material sphere_mat( (WHITE-RED)*.2, (WHITE-RED)*.10, (WHITE-RED)*.1, 12.0);
material ground_mat( GRAY*.0, GRAY*.1, GRAY*.2, 1.0);

light* lights[NUM_LIGHTS] = {
	new light(vec(1,1,1), vec(.5,.5,.5), vec(0,5,-5), 1),
	new light(vec(20,20,20), vec(5,5,5), vec(-0,5,0), 0),
};

scene_object* objs[NUM_OBJECTS] = {
	new mandlebulb(vec(1.25,0,3), sphere_mat,8,vec(0,75,75)),
	new mandlebulb(vec(-1.25,0,3), sphere_mat,6,vec(75,75,0)), 
};

vec AMBIENT(.5, 1, 1)  ;
vec SKY(75,75,50);

//TODO replace with mat defining different layers
vec pixel;


bool march(vec direction, vec ray, bool secondary = false){
	int step = 0;
	while(step < MAX_NUM_STEPS && ray.mag() < MAX_DIST){
		int index = 0;
		double min_dist = MAX_STEP;
		for( int i = 0; i < NUM_OBJECTS; i ++ ) { 
			if(objs[i]->max_dist(ray) < min_dist){
				index = i;
				min_dist = objs[i]->max_dist(ray);
			}
		}
		if( min_dist < MIN_DIST ) {
			if(!secondary) {pixel = BLACK; pixel = pixel + objs[index]->ka*AMBIENT; }
			for ( int l = 0; l < NUM_LIGHTS; l ++){
				double falloff = 1;
				vec light_direction = lights[l]->loc.unit();
				if(lights[l]->type == 0) {
					falloff = pow((lights[l]->loc-ray).mag(), -2)*10;
					light_direction = ((lights[l]->loc)-ray).unit();
				}
				vec viewer_bounce = ray.unit()-(objs[index]->normal(ray)*2*dot(ray.unit(), objs[index]->normal(ray)));
				vec begin = ray+(objs[index]->normal(ray).unit()*SHADOW_CORRECTION);
				vec shadow_direction = light_direction;
				bool render = true;
				/*if(!secondary){ //comment out this line get recursive shadows. it's really cool!!!
				if(march(shadow_direction, begin, true)) { 
					pixel = BLACK;
					pixel = pixel + objs[index]->ka*AMBIENT;
					render = false; 
				}
				}*/
				/*
				if(render && !secondary){
					pixel = pixel + objs[index]->kd*lights[l]->id*std::max(dot(light_direction, objs[index]->normal(ray)), 0.0)*falloff;
					pixel = pixel + objs[index]->ks*lights[l]->is*std::pow(std::max(dot(light_direction, viewer_bounce), 0.0), objs[index]->a)*falloff;
				}
				*/
				pixel = objs[index]->ka;
			}
			return true;
		}

		ray.step(direction,min_dist);
		step ++;
	}
	
	return false;
}


int main(){

	std::ofstream out("mandlebulb.ppm");
	out << "P3\n" << W << '\n' << H << '\n' << "255\n";
	double num = 0;
	for(double j = 0; j < H; j++){
		for(double i = 0; i < W; i ++) {
			vec ray((i-W/2)*(FRUSTUM_WIDTH/W), -(j-H/2)*(FRUSTUM_HEIGHT/H), FOCAL_LENGTH);
			pixel = BLACK;
			if( !march(ray, ORIGIN)) {
				pixel = SKY;
			}
			out << get_rgb(pixel.x) << std::endl;
			out << get_rgb(pixel.y) << std::endl;
			out << get_rgb(pixel.z) << std::endl;
			num++;
			std::cout << num/(W*H) << std::endl;
		}
	}
	std::ofstream close();
	return 0;
}
