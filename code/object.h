#ifndef OBJ_H_
#define OBJ_H_

#include"global.h"
#include"material.h"
struct Ray
{
	Point3d start;
	Vec3d direction;
	Ray() = default;
	Ray(Point3d _start, Point3d _end)
	{
		start = _start;
		direction = normalize(_end - _start);
	}
};

struct HitInfo;

class Object
{
public:
	Material material;
	Object(Material&);
	virtual void rayHit(Ray&, HitInfo&) = 0;
	virtual void computeTextureCoor(const Point3d&, double& u, double& v) = 0;
};

class Sphere :public Object
{
	double radius;
	Point3d center;
	Point3d top;
	Vec3d computeBumpNormal(const Point3d& hitPoint);
public:
	Sphere(Point3d&, double, Material&);
	void rayHit(Ray&, HitInfo&);
	void computeTextureCoor(const Point3d&, double& u, double& v);
};

class Plane :public Object			//µ¥Ãæ
{
	Point3d point;
	Vec3d axesU, axesV;  //for texture coordinates
	Vec3d normal;
public:
	Plane(Point3d&, Vec3d&, Material&, Vec3d& =Vec3d(), Vec3d& =Vec3d());
	void rayHit(Ray&, HitInfo&);
	void computeTextureCoor(const Point3d&, double& u, double& v);
};

class Rect :public Object
{
public:
	Vec3d axes[3];
	double lx, ly;
	Point3d origin;
	Rect(Point3d& o, Vec3d& vecx, Vec3d& vecy, double x, double y, Material& mat);
	void rayHit(Ray&, HitInfo&);
	void computeTextureCoor(const Point3d&, double& u, double& v);
};

class CubicBezierWall :public Object		
{
	vector<Point2d> controlPoints;
	vector<Point2d> equationCoeff;
	Object* plane;
	double height;
	void solveCubicEquation(double a, double b, double c, double d, double *ans);
	double CubicBezierWall::intersect(const Point2d& e, const Point2d& d, Point2d* intersection);
	void findPlane();
public:
	CubicBezierWall(const Point2d& p0, const Point2d& p1, const Point2d& p2, const Point2d& p3, double h, Material &mat);
	~CubicBezierWall(){ delete plane; }
	void computeTextureCoor(const Point3d&, double& u, double& v);
	void rayHit(Ray&, HitInfo&);
};

struct HitInfo
{
	double dist = INFINITY;
	Object* hitObj = NULL;
	Vec3d normal;
};

#endif