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
#include "Tool.h"

#define AllTypes 210
#define ContourNum  100
#define AngleThres 165
#define SCDisThres 0.35
#define WindowsWidth 10
#define DefaultPath "C:/Users/liuxk/OneDrive/Recent/DualNPR/"
#define SavePath "D:/vs2015project/Dihedral_new/result/"

using namespace cv;
using namespace std;

namespace Tiling_tiles {

	enum { general_p, cand_p, fea_p, fixed_p};

	typedef class Point_feature
	{
	public:
		Point2f point;
		int type; //0:普通点 1:候选点 2:特征点 3:固定点
	
	public:
		Point_feature() {};
		Point_feature(Point2f p, int t)
		{
			point = p;
			type = t;
		}
	}Point_f;

	inline Point_f p2fea(Point2f p, int type)
	{
		Point_f pf(p, type);
		return pf;
	}
	inline Point2f fea2p(Point_f p) { return p.point; }

	inline vector<Point2f> conf_trans(vector<Point_f> con)
	{
		vector<Point2f> c1;
		FOR(i, 0, con.size())
		{
			c1.push_back(con[i].point);
		}
		return c1;
	}
	typedef class innerPat
	{
	public:
		innerPat() {};
		innerPat(vector<Point_f> a, vector<int> b, int c) : in_contour(a), in_interval(b), type(c)	{};
		~innerPat() {};

		int type; //0:trans,1:rota,2:flip(13),3:flip(24)
		vector<int> in_interval;
		vector<Point_f> in_contour;
	}inPat;

	class protoTile {

	public:
		protoTile() {};
		protoTile(string filepath);
		protoTile(vector<Point2f> con);
	
		//void set_contour(vector<Point2f> c); // 给定轮廓并重采样
		//void set_anchors(vector<int> anchor_p);
		//void show_contour(vector<Point2f> c, vector<int> anchor_p);
		void feature(int  n_min, int n_max, double angle_cos);
		void resample(int sam_num);
		void scale(double scale);
		vector<Point_f>  set_flags(vector<Point2f> con, vector<int> fea);

		void Trans_contour(vector<Point_f> &c1, Point2f trans_shift);
		void Rotate_contour(vector<Point_f>& c1, Point2f center, double angle);
		void Flip_contour(vector<Point_f>& c1);
		vector<int> getCandPoints(vector<Point_f> &contour_f);


		vector<Point2f> contour;
		vector<Point_f> contour_f;
		vector<int> anchor_points;
		vector<int> feature_points;
		vector<vector<Point2f>> textures;
	};


	class Tiling_opt {

	public:
		Tiling_opt() {};
		~Tiling_opt() {};
	
		void tiliing_generation(string nameid);
		void load_dataset(bool input_images);
		
		//tiling rules
		int Tanslation_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname);
		int Flipping_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname);
		bool translation_placement(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname);

	public:
		//double dis[202][202];//两组点之间的坐标差异
		//double dis_cur[202][202];//两组点之间的曲率差异
		//double distance[202][202];
		//int step[202][202];//记录总的步数

		int all_types;
		int sampling_num;
		int allnum_inner_c;
		int match_window_width;
		int tolerance;
		double morph_ratio;
		protoTile prototile_first;
		protoTile prototile_mid;
		protoTile prototile_second;
		protoTile prototile_tem;

		vector<vector<Point_f>> contour_dataset;
		vector<vector<vector<double>>> all_con_tars;
		vector<vector<vector<double>>> all_con_tars_flip;
		vector<double> all_shape_complexity;
		vector<inPat> all_inner_conts;
		//vector<int> candadite_indexs;
		vector<vector<Point2f>> candidate_contours;
		vector<vector<pair<int, int>>> cand_paths;
	};





}