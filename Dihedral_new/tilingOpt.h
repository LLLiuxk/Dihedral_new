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

#define ClockWise 1
#define AntiClockWise 0
#define AllTypes 210
#define ContourNum  100
#define AngleThres 165
#define SCDisThres 0.35
#define WindowsWidth 10
#define DefaultPath "C:/Users/liuxk/OneDrive/Recent/DualNPR/"
#define SavePath "D:/vs2015project/Dihedral_new/result/"
#define ParaPath "D:/vs2015project/Dihedral_new/Dihedral_new/"

extern string image_id;


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
	inline vector<Point_f>  set_flags(vector<Point2f> con, vector<Point_f> con_ori)
	{
		vector<Point_f> conf;
		int csize = con.size();
		if (csize != con_ori.size())
		{
			cout << "Cannot set for con with different size!" << endl;
		}
		else
		{
			FOR(i, 0, csize)
			{
				conf.push_back(p2fea(con[i], con_ori[i].type));
			}
		}
		return conf;
	}

	typedef class innerPat
	{
	public:
		innerPat() {};
		innerPat(vector<Point_f> a, vector<int> b, int c) : in_contour(a), in_interval(b), type(c)	{
			FOR(i, 0, 4) in_contour[in_interval[i]].type = fixed_p;
		};
		~innerPat() {};

		int type; //0:trans,1:rota,2:flip(13),3:flip(24)
		vector<int> in_interval;
		vector<Point_f> in_contour;
	}inPat;

	class protoTile {

	public:
		protoTile() {};
		protoTile(string filepath, bool show=true);
		protoTile(vector<Point2f> con);
		protoTile(vector<Point_f> conf);
	
		//void set_contour(vector<Point2f> c); // 给定轮廓并重采样
		//void set_anchors(vector<int> anchor_p);
		//void show_contour(vector<Point2f> c, vector<int> anchor_p);
		void feature(int  n_min, int n_max, double angle_cos, bool show = true);
		void resample(int sam_num, bool show = true);
		void scale(double scale);
		vector<Point_f>  set_flags(vector<Point2f> con, vector<int> fea);

		void Trans_contour(vector<Point_f> &c1, Point2f trans_shift);
		void Rotate_contour(vector<Point_f>& c1, Point2f center, double angle);
		void Flip_contour(vector<Point_f>& c1);
		vector<int> getCandPoints(vector<Point_f> &contour_f);
		vector<int> getFeatures(vector<Point_f> &contour_f);
		void show_contour(vector<Point2f> c, vector<int> anchor_p);

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
		void tiliing_gen_specify(string nameid, vector<int> anc_points);
		void load_dataset(bool input_images);
		void load_para(string filename);
		void RotationVis(protoTile c1, protoTile c2, int clockorder);

		//tiling rules
		int Tanslation_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname);
		int Tanslation_rule_spec(vector<int> cand_points, vector<Point_f> &contour_s, string rootname, vector<int> anc_points);
		int Flipping_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname);
		bool translation_placement(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname);
		bool translation_spec(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname);

		void match_candidate(int inner_index);
		void match_candidate_(int inner_index);

		vector<pair<int, bool>> compare_TAR(vector<Point_f> contour_mid, int chosen_num, int window_width = WindowsWidth); //chosen_num  选择的最终结果的数目
		void feature_match(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<int> first_fea, vector<int> sec_fea, vector<pair<int, int>>& path, int &sec_shift, int width);
		vector<Point_f> morphing(vector<Point_f> contour1, vector<Point_f> contour2, vector<pair<int, int>> final_pair, double ratio);
		vector<Point_f> morph_segment(vector<Point_f> seg1, vector<Point_f> seg2, Point_f start, double ratio, double &num_error); //start 是上一个的尾端，是固定的
		vector<Point_f> morphing_dir(vector<Point_f> c_mid, vector<Point_f> c_cand, vector<pair<int, int>> path, double ratio);

	public:
		//double dis[202][202];//两组点之间的坐标差异
		//double dis_cur[202][202];//两组点之间的曲率差异
		//double distance[202][202];
		//int step[202][202];//记录总的步数

		int all_types;
		int sampling_num;
		int match_window_width;
		int tolerance;
		int Clock_order;
		double morph_ratio;
		protoTile prototile_first;
		protoTile prototile_mid;
		protoTile prototile_second;
		protoTile prototile_fin;

		vector<vector<Point_f>> contour_dataset;
		vector<vector<vector<double>>> all_con_tars;
		vector<vector<vector<double>>> all_con_tars_flip;
		vector<double> all_shape_complexity;
		vector<inPat> all_inner_conts;
		//vector<int> candadite_indexs;
		vector<vector<Point_f>> candidate_contours;
		vector<vector<pair<int, int>>> cand_paths;
		vector<vector<pair<int, int>>> cand_fea_paths;
		vector<int> path_shift;
	};





}