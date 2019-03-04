#ifndef LIGHT_H_
#define LIGHT_H_

#include"global.h"
#include"object.h"

class Light
{
public:
	Vec3d power;
	vector<Point3d> samplePoints;
	Light(Vec3d _power) :power(_power)
	{

	}
	virtual Object* object() = 0;
	virtual Ray generateRay() = 0;
};

class RectAreaLight :public Light
{
	Rect area;	
public:
	RectAreaLight(double _lx, double _ly, Point3d _origin, Vec3d vecx, Vec3d vecy, Vec3d _power) : area(_origin,vecx,vecy,_lx,_ly,Material()),Light(_power)
	{
		Vec3d p;
		Vec3d vx = area.lx / 10 * area.axes[0];
		Vec3d vy = area.ly / 10* area.axes[1];
		for (int i = 0; i < 10; i++)
		{
			p = area.origin + vx*i;
			for (int j = 0; j < 10; j++, p += vy)
			{
				samplePoints.push_back(p);
			}
		}
	}
	virtual Ray generateRay();
	virtual Object* object();
};

#endif