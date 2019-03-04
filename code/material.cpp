#include"material.h"


Material::Material(MatType _type, Vec3d _color, double _refr) :type(_type), color(_color), refractRatio(_refr)
{
	R0 = (refractRatio - 1)*(refractRatio - 1) / ((refractRatio + 1)*(refractRatio + 1));
}
Material::Material(string imgfile) : type(Diffuse)
{
	color = Vec3d(0.9, 0.9, 0.9);
	Mat img = cv::imread(imgfile);
	Vec3f tv;
	if (img.rows > 0)
	{
		size = std::min(img.rows, img.cols);
		texture = new Mat(size, size, CV_32FC3);
		Vec3f* ptr1 = &texture->at<Vec3f>(0, 0);
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				tv = (Vec3f)(img.at<Vec3b>(i, j));
				(*ptr1).val[0] = tv.val[2]/255;
				(*ptr1).val[1] = tv.val[1]/255;
				(*ptr1).val[2] = tv.val[0]/255;
				ptr1++;
			}
		}
	}
}

Material::Material(string imgfile, Vec3d _color) :type(Diffuse), color(_color)
{
	Mat bump = cv::imread(imgfile,0);
	size = std::min(bump.rows, bump.cols);
	normalMap = new Mat(size, size, CV_64FC3);
	Vec3d *ptr3 = &normalMap->at<Vec3d>(0, 0);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (j != 0 && j != size - 1)
			{
				ptr3->val[0] = -(bump.at<uchar>(i, j - 1) - bump.at<uchar>(i, j + 1))/bumpScale1;
			}
			else if (j == 0)
			{
				ptr3->val[0] = -(bump.at<uchar>(i, j) - bump.at<uchar>(i, j + 1))/bumpScale2;
			}
			else
			{
				ptr3->val[0] = -(bump.at<uchar>(i, j - 1) - bump.at<uchar>(i, j)) / bumpScale2;
			}
			if (i != 0 && i != size - 1)
			{
				ptr3->val[1] = -(bump.at<uchar>(i - 1, j) - bump.at<uchar>(i + 1, j))/bumpScale1;
			}
			else if (i == 0)
			{
				ptr3->val[1] = -(bump.at<uchar>(i, j) - bump.at<uchar>(i + 1, j)) / bumpScale2;
			}
			else
			{
				ptr3->val[1] = -(bump.at<uchar>(i - 1, j) - bump.at<uchar>(i, j))/bumpScale2;
			}
			ptr3->val[2] = 1;
			normalize(*ptr3);
			ptr3++;
		}
	}
}


Vec3d Material::getColor(double u, double v)
{
	return texture->at<Vec3f>(u*size-0.5, v*size-0.5);
}

Vec3d Material::getNormal(double u, double v)
{
	return normalMap->at<Vec3d>(u*size-0.5, v*size-0.5);
}