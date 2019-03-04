#ifndef MAP_H_
#define MAP_H_

#include"global.h"

struct Photon
{
	Point3d pos;
	Vec3d direction;
	Vec3d power;
	short layer;
	double distToPoi_2;
	
	Photon() = default;
	Photon(Point3d& _pos, Vec3d& _direction, Vec3d &_power) :pos(_pos), direction(_direction), power(_power)
	{

	}
};

class HeapCompare
{
public:
	HeapCompare()
	{

	}
	bool operator()(Photon *a, Photon *b)
	{
		return a->distToPoi_2 < b->distToPoi_2;
	}
};

class PhotonMap
{
	bool(*compare[3])(Photon&, Photon&);
	int size;
	int estimateNum;
	Point3d poi;
	double searchDist_2;

	void construction(vector<Photon>::iterator beginIt, vector<Photon>::iterator endIt, int layer, int pos);
public:
	vector<Photon> tree;
	PhotonMap(vector<Photon>& data,int estimate);
	vector<Photon*> neighbours;
	double findNeighbours(Point3d p, double dist_2);
	void rangeSearch(int pos);
};

bool compareX(Photon &a, Photon &b);
bool compareY(Photon &a, Photon &b);
bool compareZ(Photon &a, Photon &b);

#endif