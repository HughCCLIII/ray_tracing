#ifndef MATERIAL_H_
#define MATERIAL_H_
#include"global.h"

enum MatType{ Diffuse, Reflect, Refract };

struct Material
{
	Vec3d color;
	double refractRatio;
	double R0;
	MatType type;
	Mat *texture = NULL;
	Mat *normalMap = NULL;
	int size;
	Material() = default;
	Material(MatType, Vec3d, double=0);
	Material(string);
	Material(string, Vec3d);
	Vec3d getColor(double u, double v);
	Vec3d Material::getNormal(double u, double v);
};



#endif