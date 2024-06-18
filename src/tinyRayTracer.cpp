#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"
#include <numbers>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

const int width = 1024;
const int height = 768;
struct Light
{
	Light(const Vec3f& p, const float& i) : postion(p), intensity(i) {};

	Vec3f postion;
	float intensity;
};
struct Material
{
	Material(const Vec3f& color) :diffuse_color(color) {};
	Material():diffuse_color() {};
	Vec3f diffuse_color;

};
struct Sphere
{
	Vec3f center;
	float radius;
	Material material;
	Sphere(const Vec3f& c, const float& r, const Material& m) : center(c), radius(r), material(m) {}

	bool ray_intersect(const Vec3f& o, Vec3f& dir, float& t ) {
		Vec3f orig = o - center;
		double a = std::pow(dir[0], 2) + std::pow(dir[1], 2) + std::pow(dir[2], 2);
		double b = 2 * (orig[0] * dir[0] + orig[1] * dir[1] + orig[2] * dir[2]);
		double c = std::pow(orig[0], 2) + std::pow(orig[1], 2) + std::pow(orig[2], 2) - std::pow(radius, 2);
		double discr = std::pow(b, 2) - 4 * a * c;

		if (discr < 0) return false;//no intersection

		t = ((-b) - std::sqrt(discr)) / (2 * a);
		double t2 = ((-b) + std::sqrt(discr)) / (2 * a);

		if (t < 0) t = t2;//checking to see if at least one is positive
		if (t < 0) return false;

		return true;

	}

};

bool sceneIntersect(const Vec3f& orig, Vec3f& dir, std::vector<Sphere>& spheres, Vec3f& hit, Vec3f& N, Material& mat)
{
	float prv_sphere_dist = std::numeric_limits<float>::max();

	for (auto s : spheres)//for every sphere get the smallest pos tVal
	{
		float tVal=0;//scalar value for the distance to instersect

		if (s.ray_intersect(orig, dir, tVal) && tVal < prv_sphere_dist)//if tVal is smaller
		{																//(determined by raIntersect()), then that is the pixel that should be in front
			prv_sphere_dist = tVal;

			hit = orig + dir * tVal;//calculate the intersect pos//used for lighting
			N = (hit - s.center).normalize();//re orients the hit pos and noramlizes the vector
			mat = s.material;//sets material of the pixel to hit
		}
	}
	
	if (prv_sphere_dist == std::numeric_limits<float>::max())//if no sphere pixel was detected
		return false;
	
		return true;//else, prvSphereDist is different than max so there was an intersect
		

}

Vec3f cast_ray(const Vec3f& orig, Vec3f& dir, std::vector<Sphere>& spheres, std::vector<Light>& lights)
{
	Vec3f hit, N;//hit pos and normalized center
	Material mat;//material of the pixel

	if (sceneIntersect(orig, dir, spheres, hit, N, mat))//check if intersect
	{
		float diffuseLightIntensity = 0;

		for (auto l : lights)
		{
			Vec3f lightDir = (l.postion - hit).normalize();
			diffuseLightIntensity += l.intensity * std::max(0.f, lightDir * N);
		}
		return mat.diffuse_color * diffuseLightIntensity;
	}
	else return Vec3f(0.2, 0.7, 0.8);

}


void render(std::vector<Sphere>& spheres, std::vector<Light>& lights) {

	
	std::vector<Vec3f> framebuffer(width * height);//vector of Vec3f with a red blue and green value for each pixel in the image
	float fov = std::numbers::pi_v<float> / 2;//field of view

#pragma omp parallel for
	for (size_t j = 0; j < height; j++) {
		for (size_t i = 0; i < width; i++) {
			float x = (2 * (i + 0.5) / (float)width - 1)/*maps values between -1 and 1 */ * tan(fov / 2.)/*maps values between -tan and tan(fov/2)*/ * width / (float)height;//aspect ratio
			float y = -(2 * (j + 0.5) / (float)height - 1) * tan(fov / 2.);//negative because y is inverted in the image
			
			Vec3f dir = Vec3f(x, y, -1).normalize();//direction shooting rays into, so can be unsigned since normalized
			
			framebuffer[i + j * width] = cast_ray(Vec3f (0, 0, 0)/*camera pos, negative forward, pos back*/ , dir, spheres, lights);//cast a ray for every pixel
		}
	}

	std::vector<unsigned char> imageData(width * height * 3);//vector of unsigned char with an Red blue and green value in sequence for each pixel in the image
	for (size_t i = 0; i < width * height; ++i)
	{
		imageData[i * 3] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][0])));//transfering framebuffer vector to something stbi can take in
		imageData[i * 3 + 1] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][1])));
		imageData[i * 3 + 2] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][2])));

	}
	stbi_write_jpg("out.jpg", width, height, 3, imageData.data(), 100);// takes in imageData and creates a jpg file

	
}



int main() {

	Material ivory(Vec3f(0.4, 0.4, 0.3));
	Material red_rubber(Vec3f(.3, .1, .1));

	std::vector<Sphere> spheres;
	spheres.push_back(Sphere(Vec3f(-3, 0, -16), 2, ivory));
	spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
	spheres.push_back(Sphere(Vec3f(1.5, -0.5, -18), 3, red_rubber));
	spheres.push_back(Sphere(Vec3f(7, 5, -18), 4, ivory));

	std::vector<Light>  lights;
	lights.push_back(Light(Vec3f(-20, 20, 20), 1.5));

	render(spheres,lights);
	//render(s2);
	return 0;
}
