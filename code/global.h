#ifndef GLOBAL_H_
#define GLOBAL_H_

#define Debug

#include<opencv2\core.hpp>
#include<opencv2\highgui.hpp>
#include<vector>
#include<cstdio>
#include<iostream>
#include<cmath>
#include<algorithm>
#include<queue>
#include<complex>
#include<string>
#include<cassert>

using std::complex;
using std::string;
using std::vector;
using std::sqrt;
using cv::Mat;
using cv::Vec3i;
using cv::Vec3f;
using cv::Vec3b;
using cv::Vec4d;
using cv::Vec4i;
using cv::Vec2d;
using cv::Vec3d;

typedef Vec3d Point3d;
typedef Vec2d Point2d;

#define episilon 0.00001
#define PhotonNum 500000
#define CausticPhotonNum 500000
#define EstimateNum  2000
#define CausticEstimateNum 500
#define Width 1200
#define Height 1200
#define PhotonTraceDepth 6
#define RenderDepth 5
#define PhotonMapSearchDist 10
#define CausticSearchDist 10
#define FilterConstant 5.0
#define Uniform01 (rand()*1.0/RAND_MAX)
#define Uniform1_1 (Uniform01*(2*(rand()%2)-1))
#define PI 3.14159265358979
#define PhongPower 20
#define TextureScale 100
#define RayTracingSampleNum 20
const double bumpScale1 = 1.5;
const double bumpScale2 = bumpScale1 / 2;
const double filterNorm = 1 - 2.0 / (3 * FilterConstant);
const double gammaPow = 1 / 2.0;


Vec3d normalize(Vec3d& vec);

double norm(Vec3d &vec);

Vec2d normalize(Vec2d& vec);

double norm(Vec2d &vec);

int log2floor(int);

double clamp(double x); 

int gamma(double x); 

double powerWrap(double a, double k);

template<typename T>
T mysquare(T a)
{
	return a*a;
}

#endif