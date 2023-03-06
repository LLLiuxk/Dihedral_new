#pragma once

#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <ctime>
#include <io.h>
#include <list>
#include <direct.h>

using namespace cv;
using namespace std;

namespace Tiling_tiles {

	typedef struct Point_feature
	{
		Point2f point;
		int type; //0:普通点 1:候选点 2:特征点 3:固定点
	}Point_f;

	class Prototile {

	public:
		Prototile();
		Prototile(string rname, string tpath);
	
	};


	class Tiling_opt {

	public:
		Tiling_opt();
		~Tiling_opt();
	
	public:
		double dis[202][202];//两组点之间的坐标差异
		double dis_cur[202][202];//两组点之间的曲率差异
		double distance[202][202];
		int step[202][202];//记录总的步数

		int all_types;
		int sampling_num;
		int allnum_inner_c;
		int match_window_width;
		int tolerance;
		double morph_ratio;
		Prototile *prototile_first;
		Prototile *prototile_mid;
		Prototile *prototile_second;
		Prototile *prototile_tem;
	};





}