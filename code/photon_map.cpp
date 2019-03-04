#include"photon_map.h"

bool compareX(Photon &a, Photon &b)
{
	return a.pos.val[0] < b.pos.val[0];
}

bool compareY(Photon &a, Photon &b)
{
	return a.pos.val[1] < b.pos.val[1];
}

bool compareZ(Photon &a, Photon &b)
{
	return a.pos.val[2] < b.pos.val[2];
}


PhotonMap::PhotonMap(vector<Photon>& data, int estimate) :estimateNum(estimate)
{
	compare[0] = compareX;
	compare[1] = compareY;
	compare[2] = compareZ;
	size = data.size() + 1;
	tree.resize(size);
	construction(data.begin(),data.end(), 0, 1);
	if (tree.size() % 2 == 1)
	{
		Photon p(Vec3d(0, 0, 0), Vec3d(0, 0, 0), Vec3d(0, 0, 0));
		tree.push_back(p);
	}
	size = tree.size();
	tree.shrink_to_fit();
}

void PhotonMap::construction(vector<Photon>::iterator beginIt, vector<Photon>::iterator endIt, int layer, int pos)
{
	int treeSize = endIt - beginIt + 1;
	int tsize = 1 << (log2floor(treeSize));
	tsize = treeSize >= tsize + (tsize >> 1) ? tsize : (tsize >> 1);
	int leftTreeSize = std::max(tsize - 1, treeSize - tsize - 1);
	if (endIt - beginIt == 1)
	{
		tree[pos] = *beginIt;
		tree[pos].layer = layer;
		return;
	}
	nth_element(beginIt, beginIt + leftTreeSize, endIt,compare[layer]);
	tree[pos] = *(beginIt + leftTreeSize);
	tree[pos].layer = layer;
	layer = (layer + 1) % 3;
	auto partitionIt = beginIt + leftTreeSize;
	construction(beginIt, partitionIt, layer, 2 * pos);
	if (endIt > partitionIt + 1)
	{
		construction(partitionIt + 1, endIt, layer, 2 * pos + 1);
	}	
}

double PhotonMap::findNeighbours(Point3d p, double dist_2)
{
	neighbours.clear();
	poi = p;
	searchDist_2 = dist_2;
	rangeSearch(1);
	if (!neighbours.empty())
		return (*neighbours.begin())->distToPoi_2;
	return 1;
}


void PhotonMap::rangeSearch(int pos)
{
	Photon *p = &tree[pos];
	if (2 * pos + 1 <= size)
	{
		int layer = p->layer;
		double d = poi.val[layer] - p->pos.val[layer];
		if (d < 0)
		{
			rangeSearch(2 * pos);
			if (d*d < searchDist_2)
			{
				rangeSearch(2 * pos + 1);
			}
		}
		else
		{
			rangeSearch(2 * pos + 1);
			if (d*d < searchDist_2)
			{
				rangeSearch(2 * pos);
			}
		}
	}
	Vec3d t = poi - p->pos;
	p->distToPoi_2 = t.ddot(t);
	if (p->distToPoi_2 < searchDist_2)
	{
		neighbours.push_back(p);
		std::push_heap(neighbours.begin(), neighbours.end(), HeapCompare());
		if (neighbours.size() > estimateNum)
		{
			std::pop_heap(neighbours.begin(), neighbours.end(), HeapCompare());
			neighbours.pop_back();
			searchDist_2 = neighbours[0]->distToPoi_2;
		}
	}	
}



