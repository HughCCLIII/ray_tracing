#include"scene.h"

Camera::Camera(Point3d position, Point3d focus, double dist) :pos(position), imageDist(dist)
{
	axes[1] = normalize(focus - pos);
	axes[0] = normalize(axes[1].cross(Vec3d(0, 0, 1)));
	axes[2] = axes[0].cross(axes[1]);
}

Scene::Scene(Camera _cam, Light* _light) :cam(_cam), light(_light)
{

}

void Scene::addObj(Object* obj)
{
	objs.push_back(obj);
}

void Scene::findFirstObjectHit(Ray& ray, HitInfo& hitInfo)
{
	for (auto it = objs.begin(); it != objs.end(); it++)
	{
		(*it)->rayHit(ray, hitInfo);
	}
}

void Scene::render(Mat& img)
{
	double gridSize = std::min(2.0 / Width, 2.0 / Height);
	Point3d leftUp(-gridSize*(Width - 1) / 2, cam.imageDist, gridSize*(Height - 1) / 2);
	Point3d canvasLeftUp = leftUp.val[0] * cam.axes[0] + leftUp.val[1] * cam.axes[1] + leftUp.val[2] * cam.axes[2] + cam.pos;
	Vec3d gridX = gridSize * cam.axes[0];
	Vec3d gridZ = gridSize * cam.axes[2];
	Ray ray;
	Material m;
	int depth;
	Vec3d color;
	Vec3d sampleColor;
	Point3d hitPoint;
	Point3d samplePoint;
	Vec3d halfgridX = gridX / 2;
	Vec3d halfgridZ = gridZ / 2;
	for (int i = 0; i < Height; i++)
	{
		Point3d renderPixel = canvasLeftUp - i * gridZ;
		for (int j = 0; j < Width; j++, renderPixel+=gridX)
		{
			color.val[0] = color.val[1] = color.val[2] = 0;
			for (int k = 0; k < RayTracingSampleNum; k++)
			{
				HitInfo hitInfo;
				sampleColor.val[0] = sampleColor.val[1] = sampleColor.val[2] = 0;
				samplePoint = renderPixel + Uniform1_1*halfgridX + Uniform1_1*halfgridZ;
				ray.start = cam.pos;
				ray.direction = normalize(samplePoint - cam.pos);
				depth = 0;

				findFirstObjectHit(ray, hitInfo);
				if (hitInfo.hitObj)
				{
					sampleColor = computeColor(ray.start + hitInfo.dist * ray.direction, ray.direction, hitInfo, 0);
				}
				color += sampleColor / RayTracingSampleNum;
			}
			img.at<Vec3b>(i, j) = Vec3b(gamma(color.val[2]), gamma(color.val[1]), gamma(color.val[0]));
		}
	}
	
}

void Scene::generatePhotonMap()
{
	for (int i = 0; i < PhotonNum; i++)
	{
		while (!emitPhoton(light));
	}
	globalMap = new PhotonMap(globalPhotons, EstimateNum);
	globalPhotons.clear();
	globalPhotons.shrink_to_fit();
	while (causticPhotons.size() < CausticPhotonNum)
	{
		emitCausticPhoton(light);
	}
	causticMap = new PhotonMap(causticPhotons, CausticEstimateNum);
	causticPhotons.clear();
	causticPhotons.shrink_to_fit();
}

void Scene::emitCausticPhoton(Light *light)
{
	int depth = 0;
	Vec3d power = light->power;
	Ray ray;
	ray = light->generateRay();
	Material m;
	Photon photon;
	while (depth < PhotonTraceDepth)
	{
		HitInfo hitInfo;
		findFirstObjectHit(ray, hitInfo);
		if (!hitInfo.hitObj || hitInfo.hitObj == light->object())
			return;
		m = hitInfo.hitObj->material;
		switch (m.type)
		{
		case Diffuse:
			if (depth)
			{
				ray.start += hitInfo.dist*ray.direction;
				power = power.mul(m.color);
				photon.power = power;
				photon.direction = ray.direction;
				photon.pos = ray.start;
				causticPhotons.push_back(photon);
				return;
			}
			else
				return;
			break;
		case Reflect:
			power = power.mul(m.color);
			ray.start += hitInfo.dist*ray.direction;
			computeReflectDirection(hitInfo.normal, ray.direction, ray.direction);
			break;
		case Refract:
		{
			bool inside;
			ray.start += (hitInfo.dist)*ray.direction;		//保证点进入物体内部
			if (checkAllReflection(hitInfo.normal, ray.direction, m.refractRatio, inside))
			{
				computeReflectDirection(-hitInfo.normal, ray.direction, ray.direction);
				power = power.mul(m.color);
			}
			else
			{
				double randp = Uniform01;
				double reflectance, c;
				if (inside)
				{
					c = sqrt(1 - mysquare<double>(m.refractRatio)*(1 - mysquare<double>(hitInfo.normal.ddot(ray.direction))));
					c = 1 - c;
				}
				else
				{
					c = 1 + hitInfo.normal.ddot(ray.direction);
				}
				reflectance = m.R0 + (1 - m.R0)*c*mysquare<double>(mysquare<double>(c));
				if (randp < reflectance)
				{
					computeReflectDirection((inside?-1:1)*hitInfo.normal, ray.direction, ray.direction);
					power = power.mul(m.color);
				}
				else
				{
					computeRefractDirection(hitInfo.normal, ray.direction, ray.direction, m.refractRatio);
					power = power.mul(m.color);
				}				
			}
			break;
		}		
		}
		depth++;
	}
}

bool Scene::emitPhoton(Light *light)
{
	int depth = 0;
	Vec3d power = light->power;
	Ray ray;
	ray = light->generateRay();
	Material m;
	Photon photon;
	bool success = false;
	while (depth < PhotonTraceDepth)
	{
		HitInfo hitInfo;
		findFirstObjectHit(ray, hitInfo);
		if (!hitInfo.hitObj || hitInfo.hitObj == light->object())
			return success;
		m = hitInfo.hitObj->material;
		switch (m.type)
		{
		case Diffuse:
		{
			ray.start += hitInfo.dist*ray.direction;
			power = power.mul(m.color);
			photon.power = power;
			photon.direction = ray.direction;
			photon.pos = ray.start;
			globalPhotons.push_back(photon);
			success = true;

			Vec3d basis1 = normalize(hitInfo.normal.cross(Vec3d(rand(), rand(), rand())));
			Vec3d basis2 = hitInfo.normal.cross(basis1);
			double xx, yy, zz;
			do{
				xx = Uniform1_1;
				yy = Uniform1_1;
				zz = Uniform01;
			} while (xx*xx + yy*yy + zz*zz > 1);
			ray.direction = normalize(xx*basis1 + yy*basis2 + zz*hitInfo.normal);
			break;
		}

		case Reflect:
			power = power.mul(m.color);
			ray.start += hitInfo.dist*ray.direction;
			computeReflectDirection(hitInfo.normal, ray.direction, ray.direction);
			break;
		case Refract:
		{
			bool inside;
			ray.start += (hitInfo.dist)*ray.direction;		//保证点进入物体内部
			if (checkAllReflection(hitInfo.normal, ray.direction, m.refractRatio, inside))
			{
				computeReflectDirection(-hitInfo.normal, ray.direction, ray.direction);
				power = power.mul(m.color);
			}
			else
			{
				double randp = Uniform01;
				double reflectance, c;
				if (inside)
				{
					c = sqrt(1 - mysquare<double>(m.refractRatio)*(1 - mysquare<double>(hitInfo.normal.ddot(ray.direction))));
					c = 1 - c;
				}
				else
				{
					c = 1 + hitInfo.normal.ddot(ray.direction);
				}
				reflectance = m.R0 + (1 - m.R0)*c*mysquare<double>(mysquare<double>(c));
				if (randp < reflectance)
				{
					computeReflectDirection((inside ? -1 : 1)*hitInfo.normal, ray.direction, ray.direction);
					power = power.mul(m.color);
				}
				else
				{
					computeRefractDirection(hitInfo.normal, ray.direction, ray.direction, m.refractRatio);
					power = power.mul(m.color);
				}
			}
			break;
		}
		}
		depth++;
	}
	return success;
}

void computeReflectDirection(const Vec3d& normal, const Vec3d& incident, Vec3d& reflect)
{
	reflect = incident - 2 * incident.ddot(normal)*normal;
}

bool checkAllReflection(const Vec3d& normal, const Vec3d& incident, double refr, bool& inside) //true表示全反射
{
	double innerProduct = normal.ddot(incident);
	double sign = 1;
	inside = innerProduct > 0;
	return (inside && innerProduct < sqrt(1 - 1 / (refr*refr)));

}

void computeRefractDirection(const Vec3d& normal, const Vec3d& incident, Vec3d& refract,  double refr)
{
	double innerProduct = normal.ddot(incident);
	double sign = 1;
	if (innerProduct > 0)   //射向空气
	{
		refr = 1 / refr;
		sign = -1;
	}
	Vec3d x = incident - innerProduct*normal;
	refract = x / refr - sqrt(1 - (1 - innerProduct*innerProduct) / (refr*refr))*normal*sign;
}


Vec3d Scene::radianceAtDiff(const Point3d& hitPoint, const Vec3d& incidence, const HitInfo& hitInfo)
{
	Vec3d lightColor;
	Ray ray;
	ray.start = hitPoint;
	Vec3d radiance = Vec3d(0.2, 0.2, 0.2) / (double)light->samplePoints.size();
	Material& m = hitInfo.hitObj->material;
	double radius_2, radius;
	Vec3d colorTemp;


	for (auto it = light->samplePoints.begin(); it != light->samplePoints.end(); it++)
	{
		HitInfo hinfo;
		ray.direction = normalize(*it - hitPoint);
		findFirstObjectHit(ray, hinfo);
		if (hinfo.hitObj == light->object())
		{
			lightColor += (radiance * std::max(0.0, ray.direction.ddot(hitInfo.normal)));
			lightColor += (radiance * std::max(0.0, std::pow(hitInfo.normal.ddot(normalize(ray.direction - incidence)), PhongPower)));
		}
	}

	Vec3d diffuseColor;
	radius_2 = globalMap->findNeighbours(hitPoint, PhotonMapSearchDist);
	radius = sqrt(radius_2);
	vector<Photon*> &neighbours = globalMap->neighbours;
	for (auto it = neighbours.begin(); it != neighbours.end(); it++)
	{
		colorTemp += (*it)->power * std::max(0.0,(*it)->direction.ddot(-hitInfo.normal))*(1 - norm(hitPoint - (*it)->pos) / (radius*FilterConstant));
	}
	
	diffuseColor += colorTemp / (PI*radius_2*filterNorm);	
	colorTemp = Vec3d(0, 0, 0);
	radius_2 = causticMap->findNeighbours(hitPoint, CausticSearchDist);
	radius = sqrt(radius_2);
	neighbours = causticMap->neighbours;
	for (auto it = neighbours.begin(); it != neighbours.end(); it++)
	{
		colorTemp += (*it)->power * std::max(0.0, (*it)->direction.ddot(-hitInfo.normal))* (1 - norm(hitPoint - (*it)->pos) / (radius*FilterConstant));
	}
	diffuseColor += colorTemp / (PI*radius_2*filterNorm);
	

	return lightColor + diffuseColor;
}


Vec3d Scene::computeColor(const Point3d& hitPoint, const Vec3d& incidence, const HitInfo& hitInfo, int depth)
{
	if (depth >= RenderDepth)
	{
		return Vec3d(0, 0, 0);
	}
		
	if (hitInfo.hitObj == light->object())
	{
		return Vec3d(1, 1, 1);
	}
	else
	{
		Vec3d color(0,0,0);
		Material& m = hitInfo.hitObj->material;
		switch (m.type)
		{
		case Diffuse:
		{
			double u, v;
			if (m.texture)
			{
				hitInfo.hitObj->computeTextureCoor(hitPoint, u, v);
				color += m.getColor(u, v).mul(radianceAtDiff(hitPoint, incidence, hitInfo));
			}
			else
			{
				color += m.color.mul(radianceAtDiff(hitPoint, incidence, hitInfo));
			}
			break;
		}
		case Reflect:
		{
			Vec3d newDirection;
			computeReflectDirection(hitInfo.normal, incidence, newDirection);
			HitInfo hinfo;
			Ray ray;
			ray.start = hitPoint;
			ray.direction = newDirection;
			findFirstObjectHit(ray, hinfo);
			if (hinfo.hitObj)
			{
				color += m.color.mul(computeColor(hitPoint + hinfo.dist*newDirection, newDirection, hinfo, depth + 1));
			}
			break;
		}
			
		case Refract:
		{
			double reflectance;
			bool inside;
			Vec3d reflectDirection, refractDirection;
			if (checkAllReflection(hitInfo.normal, incidence, m.refractRatio, inside))
			{
				reflectance = 1;
			}
			else
			{
				double c;
				if (inside)
				{
					c = sqrt(1 - mysquare<double>(m.refractRatio)*(1 - mysquare<double>(hitInfo.normal.ddot(incidence))));
					c = 1 - c;
				}
				else
				{
					c = 1 + hitInfo.normal.ddot(incidence);
				}
				reflectance = m.R0 + (1 - m.R0)*c*mysquare<double>(mysquare<double>(c));
			}
			HitInfo hinfo1,hinfo2;
			Ray ray1,ray2;
			ray2.start = ray1.start = hitPoint;
			computeReflectDirection((inside?-1:1)*hitInfo.normal,incidence,ray1.direction);
			findFirstObjectHit(ray1, hinfo1);
			if (hinfo1.hitObj)
			{
				color += m.color.mul(computeColor(hitPoint + hinfo1.dist*ray1.direction, ray1.direction, hinfo1, depth + 1))*reflectance;
			}
			if (reflectance < 1 - episilon)
			{
				computeRefractDirection(hitInfo.normal, incidence, ray2.direction, m.refractRatio);
				findFirstObjectHit(ray2, hinfo2);
				if (hinfo2.hitObj)
				{
					color += m.color.mul(computeColor(hitPoint + hinfo2.dist*ray2.direction, ray2.direction, hinfo2, depth + 1))*(1 - reflectance);
				}
			}		
			break;
		}
			
		}
		return color;
	}
}