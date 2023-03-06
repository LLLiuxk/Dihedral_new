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
		int type; //0:��ͨ�� 1:��ѡ�� 2:������ 3:�̶���
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
		double dis[202][202];//�����֮����������
		double dis_cur[202][202];//�����֮������ʲ���
		double distance[202][202];
		int step[202][202];//��¼�ܵĲ���

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