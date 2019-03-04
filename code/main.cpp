#include "global.h"
#include"scene.h"
	
int main()
{
	Material concrete("texture/Concrete.tif");
	Material tile("texture/tile1.jpg");
	Material bump("texture/StoneWall_height.tif",Vec3d(0.8,0.8,0.8));
	Material marbleFloor("texture/MarbleFloor1.tiff");
	Material allReflect(MatType::Reflect, Vec3d(1, 1, 1));
	Material refract1(MatType::Refract, Vec3d(0.75, 0.75, 1),1.6);
	Material refract2(MatType::Refract, Vec3d(1, 0.75, 1), 1.6);
	Material refract3(MatType::Refract, Vec3d(1, 0.8, 0), 1.6);
	Material whiteWall(MatType::Diffuse, Vec3d(0.75, 0.75, 0.75));
	Camera cam(Point3d(50, -20, 50), Point3d(-100, 20, 50),1.1);
	RectAreaLight light(30, 30, Point3d(-35, -15, 100), Vec3d(-1, 0, 0), Vec3d(0, 1, 0), Vec3d(7000.0/PhotonNum, 7000.0/PhotonNum,7000.0/PhotonNum));
	Scene scene(cam,&light);
	Plane floor(Point3d(-50, 0, 0), Vec3d(0, 0, 1), marbleFloor, Vec3d(1, 0, 0), Vec3d(0, 1, 0));
	Plane ceiling(Point3d(-50, 0, 100), Vec3d(0, 0, -1), bump, Vec3d(-1, 0, 0), Vec3d(0, 1, 0));
	Plane leftWall(Point3d(-50, -70, 50), Vec3d(0, 1, 0), tile, Vec3d(-1, 0, 0), Vec3d(0, 0, 1));
	Plane rightWall(Point3d(-30, 80, 50), Vec3d(0, -1, 0), tile, Vec3d(-1, 0, 0), Vec3d(0, 0, 1));
	Plane frontWall(Point3d(10, 0, 50), Vec3d(-1, 0, 0), concrete, Vec3d(0, 1, 0), Vec3d(0, 0, 1));
	Sphere ball1(Point3d(-40, 40, 15), 15, refract1);
	Sphere ball2(Point3d(-55, -30, 15), 15, refract2);
	Sphere ball3(Point3d(-20, -2, 15), 15, refract3);
	CubicBezierWall backWall1(Point2d(-100, -70), Point2d(-75, -35), Point2d(-110, 0), Point2d(-110, 20), 100, allReflect);
	CubicBezierWall backWall2(Point2d(-110, 20), Point2d(-110, 45), Point2d(-90, 60), Point2d(-100, 80), 100, allReflect);
	scene.addObj(light.object());
	scene.addObj(&floor);
	scene.addObj(&ceiling);
	scene.addObj(&leftWall);
	scene.addObj(&rightWall);
	scene.addObj(&frontWall);
	scene.addObj(&ball3);
	scene.addObj(&backWall1);
	scene.addObj(&backWall2);
	scene.addObj(&ball1);
	scene.addObj(&ball2);
	
	Mat image(Height, Width, CV_8UC3);
	scene.generatePhotonMap();
	scene.render(image);

	cv::namedWindow("result");
	cv::imshow("result", image);
	cv::waitKey(0);
	cv::imwrite("laji14.png", image);
	return 0;

}