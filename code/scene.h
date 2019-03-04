#ifndef SCENE_H_
#define SCENE_H_

#include "global.h"
#include "photon_map.h"
#include "material.h"
#include "object.h"
#include "light.h"

struct Camera
{
	Point3d pos;
	Vec3d axes[3];
	double imageDist;
	Camera(Point3d position, Point3d focus, double dist);
};

class Scene
{
	vector<Object*> objs;
	Camera cam;
	Light *light;
	PhotonMap *globalMap, *causticMap;
	vector<Photon> globalPhotons, causticPhotons;
public:
	Scene(Camera _cam, Light *_light);

	void addObj(Object*);
	void findFirstObjectHit(Ray&, HitInfo&);
	void generatePhotonMap();
	void render(Mat&);
	bool emitPhoton(Light*);
	void emitCausticPhoton(Light *light);
	Vec3d Scene::radianceAtDiff(const Point3d& hitPoint, const Vec3d& incidence, const HitInfo& hitInfo);
	Vec3d computeColor(const Point3d& hitPoint, const Vec3d& incidence, const HitInfo& hitInfo, int depth);
};

void computeReflectDirection(const Vec3d& normal, const Vec3d& incident, Vec3d& reflect);
void computeRefractDirection(const Vec3d& normal, const Vec3d& incident, Vec3d& refract, double refr);
bool checkAllReflection(const Vec3d& normal, const Vec3d& incident, double refr, bool& inside);
#endif