#include "light.h"

Ray RectAreaLight::generateRay()
{
	Ray ray;
	ray.start = area.origin + Uniform01 * area.lx * area.axes[0]+ Uniform01 * area.ly * area.axes[1];
	double x, y, z;
	do{
		x = Uniform1_1;
		y = Uniform1_1;
		z = Uniform01;
	} while (x*x + y*y + z*z > 1);
	ray.direction = normalize(x*area.axes[0] + y*area.axes[1] + z*area.axes[2]);
	return ray;
}

Object *RectAreaLight::object()
{
	return &area;
}