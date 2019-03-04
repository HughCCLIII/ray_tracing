#include"object.h"

Object::Object(Material& mat) :material(mat)
{

}

Sphere::Sphere(Point3d& c, double r, Material& mat) :center(c), radius(r), Object(mat)
{
	top = center + Vec3d(0, 0, radius);
}

void Sphere::rayHit(Ray& ray, HitInfo& hitInfo)
{
	Vec3d tv = ray.start - center;
	double ddot1 = ray.direction.ddot(tv),ddot2=tv.ddot(tv);
	double delta = ddot1*ddot1 - ddot2 + radius*radius;
	double t1, t2;
	if (delta > 0)
	{
		delta = sqrt(delta);
		t1 = -ddot1 - delta;
		t2 = -ddot1 + delta;
		if (t1 > episilon && t1 < hitInfo.dist)
		{
			hitInfo.hitObj = this;
			hitInfo.dist = t1;
			if (material.normalMap)
			{
				hitInfo.normal = computeBumpNormal(ray.start + t1*ray.direction);
			}
			else
			{
				hitInfo.normal = normalize(ray.start + t1*ray.direction - center);
			}
			return;
		}
		if (t2 > episilon && t2 < hitInfo.dist)
		{
			hitInfo.hitObj = this;
			hitInfo.dist = t2;
			if (material.normalMap)
			{
				hitInfo.normal = computeBumpNormal(ray.start + t2*ray.direction);
			}
			else
			{
				hitInfo.normal = normalize(ray.start + t2*ray.direction - center);
			}		
		}
	}
}

Vec3d Sphere::computeBumpNormal(const Point3d& hitPoint)
{
	Vec3d axesX, axesY, axesZ;
	axesZ = normalize(hitPoint - center);
	axesX = normalize(hitPoint - top).cross(axesZ);
	axesY = axesZ.cross(axesY);
	double u, v;
	computeTextureCoor(hitPoint, u, v);
	Vec3d coor = material.getNormal(u, v);
	return coor.val[0] * axesX + coor.val[1] * axesY + coor.val[2] * axesZ;
}

void Sphere::computeTextureCoor(const Point3d& poi, double& u, double& v)
{
	Point3d tp = poi - center;
	u = 0.5 + atan2(tp.val[1], tp.val[0]) / (2 * PI);
	v = 1 - acos(tp.val[2] / radius) / PI;
}

Plane::Plane(Point3d& _p, Vec3d& _normal, Material& mat, Vec3d& _axesU, Vec3d& _axesV) : point(_p), normal(_normal), Object(mat), axesU(_axesU), axesV(_axesV)
{
	
}

void Plane::rayHit(Ray& ray, HitInfo& hitInfo)
{
	if (normal.ddot(ray.start - point) <= 0)
		return;
	double d_hit = (point-ray.start).ddot(normal) / ray.direction.ddot(normal);
	if (d_hit > episilon && d_hit < hitInfo.dist)
	{
		hitInfo.hitObj = this;
		hitInfo.dist = d_hit;
		if (material.normalMap)
		{
			double u, v;
			computeTextureCoor(ray.start + d_hit*ray.direction, u, v);
			Vec3d normalCoor = material.getNormal(u, v);
			hitInfo.normal = axesU*normalCoor.val[0] + axesV*normalCoor.val[1] + normal*normalCoor.val[2];
		}
		else
		{
			hitInfo.normal = normal;
		}
	}
}


void Plane::computeTextureCoor(const Point3d& p, double& u, double& v)
{
	Vec3d tv = p - point;
	u = tv.ddot(axesU)/TextureScale;
	v = tv.ddot(axesV)/TextureScale;
	u -= (int)u;
	u = u > 0 ? u : u + 1;
	v -= (int)v;
	v = v > 0 ? v : v + 1;
}


Rect::Rect(Point3d& o, Vec3d& vecx, Vec3d& vecy, double x, double y, Material& mat) :Object(mat), lx(x), ly(y), origin(o)
{
	axes[0] = vecx;
	axes[1] = vecy;
	axes[2] = vecx.cross(vecy);
}

void Rect::rayHit(Ray& ray, HitInfo& hitInfo)
{
	double d_hit = (origin - ray.start).ddot(axes[2]) / ray.direction.ddot(axes[2]);
	if (d_hit > episilon && d_hit < hitInfo.dist)
	{
		Point3d vecHit = ray.start + d_hit*ray.direction - origin;
		double xcoor = vecHit.ddot(axes[0]), ycoor = vecHit.ddot(axes[1]);
		if (xcoor <= lx && xcoor>=0 && ycoor <= ly && ycoor >=0)
		{
			hitInfo.hitObj = this;
			hitInfo.dist = d_hit;
			hitInfo.normal = axes[2] * (axes[2].ddot(ray.start - origin) <= 0 ? -1 : 1);
		}

	}
}

void Rect::computeTextureCoor(const Point3d&, double& u, double& v)
{

}


CubicBezierWall::CubicBezierWall(const Point2d& p0, const Point2d& p1, const Point2d& p2, const Point2d& p3, double h, Material& mat) :Object(mat), height(h)
{
	controlPoints.push_back(p0);
	controlPoints.push_back(p1);
	controlPoints.push_back(p2);
	controlPoints.push_back(p3);
	equationCoeff.push_back(3 * (p1 - p2) + p3 - p0);
	equationCoeff.push_back(3 * (p2 + p0 - 2 * p1));
	equationCoeff.push_back(3 * (p1 - p0));
	equationCoeff.push_back(p0);

	findPlane();
}

void CubicBezierWall::rayHit(Ray& ray, HitInfo& hinfo)	//先判断和包围平面是否有交 若有 投影到xy平面(墙面与xy平面垂直) 然后解方程求交
{
	HitInfo tinfo;
	plane->rayHit(ray, tinfo);
	if (tinfo.hitObj == plane && tinfo.dist < hinfo.dist)
	{
		Point2d e(ray.start.val[0], ray.start.val[1]);
		Point2d d(ray.direction.val[0], ray.direction.val[1]);
		normalize(d);
		Point2d intersection;
		double t;
		Vec2d deriv;
		t = intersect(e, d, &intersection);
		if (t > 0)
		{
			double tan = ray.direction.val[2] / sqrt(ray.direction.val[0] * ray.direction.val[0] + ray.direction.val[1] * ray.direction.val[1]);
			double hdist = norm(intersection - e);
			double h = hdist*tan + ray.start.val[2];
			if (h < height && h>0)
			{
				double dist = sqrt(hdist*hdist*(1 + tan*tan));
				if (dist < hinfo.dist)
				{
					hinfo.dist = dist;
					hinfo.hitObj = this;
					deriv = 3 * t * t * equationCoeff[0] + 2 * t*equationCoeff[1] + equationCoeff[2];
					Vec3d normal3d(-deriv.val[1], deriv.val[0], 0);
					normalize(normal3d);
					normal3d = normal3d.ddot(ray.direction) < 0 ? normal3d : -normal3d;
					hinfo.normal = normal3d;
				}
			}
		}
	}
}

double CubicBezierWall::intersect(const Point2d& e, const Point2d& d, Point2d* intersection)
{
	double ans[3];
	double t = -1;
	Vec2d dd(-d.val[1], d.val[0]);
	solveCubicEquation(dd.ddot(equationCoeff[0]), dd.ddot(equationCoeff[1]), dd.ddot(equationCoeff[2]), dd.ddot(equationCoeff[3] - e), ans);
	Point2d  tp;
	double dist = INFINITY, td;
	double c1, c2, c3;
	for (int i = 0; i < 3; i++)
	{
		if (ans[i] > 0 && ans[i] < 1)
		{
			c3 = ans[i];
			c2 = c3*c3;
			c1 = c2*c3;
			tp = equationCoeff[0] * c1 + equationCoeff[1] * c2 + equationCoeff[2] * c3 + equationCoeff[3];
			int id = 0;
			if (d.val[0] * d.val[0] < d.val[1] * d.val[1])
			{
				id = 1;
			}
			td = (tp.val[id] - e.val[id]) / d.val[id];
			if (td < dist)
			{
				dist = td;
				if(intersection)
					*intersection = tp;
				t = c3;
			}
		}
	}
	return t;
}


void CubicBezierWall::solveCubicEquation(double a, double b, double c, double d, double* ans)
{
	double p, q;
	p = (3 * a*c - b*b) / (9 * a*a);
	q = (27 * a*a*d - 9 * a*b*c + 2 * b*b*b) / (54 * a*a*a);
	double delta = q*q + p*p*p;
	double sdelta = b / (3 * a);
	if (delta > 0)
	{
		delta = sqrt(delta);
		double x1 = powerWrap(-q + delta, 0.33333) + powerWrap(-q - delta, 0.33333);

		x1 -= sdelta;
		ans[0] = x1;
		ans[1] = -1;
		ans[2] = -1;
	}
	else
	{
		complex<double> t1, t2, dd(delta);
		complex<double> o1(-0.5, sqrt(3) / 2.0),o2;
		o2._Val[0] = -0.5;
		o2._Val[1] = -o1.imag();
		dd = sqrt(dd);
		t1 = pow(-q + dd, 0.33333);
		t2 = pow(-q - dd, 0.33333);
		ans[0] = (t1 + t2).real() - sdelta;
		ans[1] = (o1*t1 + o2*t2).real() - sdelta;
		ans[2] = (o2*t1 + o1*t2).real() - sdelta;
	}
}

void CubicBezierWall::findPlane()
{
	Point2d start = controlPoints[0];
	Vec2d d;
	Vec2d p = normalize(controlPoints[3] - controlPoints[0]);
	Vec3d td = Vec3d(p.val[0], p.val[1], 0).cross(Vec3d(0, 0, 1));
	d.val[0] = td.val[0];
	d.val[1] = td.val[1];
	Point2d end = start + 50 * d;
	Point2d e;
	for (int i = 0; i < 10; i++)
	{
		e = (start + end) / 2;
		double t = intersect(e, p, NULL);
		if (t > 0 && t < 1)
		{
			start = e;
		}
		else
		{
			end = e;
		}
	}
	e += 0.1*d;
	plane = new Plane(Point3d(e.val[0], e.val[1], 0), Vec3d(d.val[0],d.val[1],0),Material());
}

void CubicBezierWall::computeTextureCoor(const Point3d&, double& u, double& v)
{

}
