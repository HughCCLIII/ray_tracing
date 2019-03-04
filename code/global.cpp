#include"global.h"

Vec3d normalize(Vec3d& vec)
{
	vec /= norm(vec);
	return vec;
}

Vec2d normalize(Vec2d& vec)
{
	vec /= norm(vec);
	return vec;
}

int log2floor(int x) {
	int bitsNumber = 0;
	bitsNumber = (!!(x >> 16)) << 4;
	bitsNumber = bitsNumber + ((!!(x >> (bitsNumber + 8))) << 3);
	bitsNumber = bitsNumber + ((!!(x >> (bitsNumber + 4))) << 2);
	bitsNumber = bitsNumber + ((!!(x >> (bitsNumber + 2))) << 1);
	bitsNumber = bitsNumber + (!!(x >> (bitsNumber + 1)));
	bitsNumber = bitsNumber + (!!bitsNumber) + (~0) + (!(1 ^ x));
	return bitsNumber;
}

double norm(Vec3d &vec)
{
	return sqrt(vec.ddot(vec));
}

double norm(Vec2d &vec)
{
	return sqrt(vec.ddot(vec));
}

double clamp(double x){ return x < 0 ? 0 : (x > 1 ? 1 : x); }

int gamma(double x) { return int(pow(clamp(x), gammaPow) * 255 + .5); }

double powerWrap(double a, double k)
{
	if (a < 0)
		return -pow(-a, k);
	return pow(a, k);
}

