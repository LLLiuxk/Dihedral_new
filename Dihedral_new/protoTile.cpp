#include "tilingOpt.h"

namespace Tiling_tiles {

	protoTile::protoTile(string filepath)
	{
		contour = load_point_file(filepath);
		if (contour.empty())
		{
			cout << filepath << " does not exist! " << endl;
		}
		else
		{
			//cout << "hahahah:" << csize << endl;
			resample(ContourNum);
			contour_f = set_flags(contour, feature_points);
		}		
	}

	protoTile::protoTile(vector<Point2f> con)
	{
		contour = con;
		resample(ContourNum);
		contour_f = set_flags(contour, feature_points);
	}


	void protoTile::feature(int  n_min, int n_max, double angle_cos)
	{
		feature_points = cal_feature(contour, n_min, n_max, angle_cos, true);
	}


	void protoTile::resample(int sam_num)
	{
		int csize = contour.size();
		feature(1, 0.02*csize, cos(PI * AngleThres / 180));
		contour = con_sample(contour, feature_points, sam_num, true);
	}


	void protoTile::scale(double scale)
	{
		Point2f cen = center_p(contour);
		int csize = contour.size();
		for (int i = 0; i < csize; i++)
		{
			contour[i] = scale * (contour[i] - cen) + cen;
		}
	}

	vector<Point_f> protoTile::set_flags(vector<Point2f> con, vector<int> fea)
	{
		vector<Point_f> conf;
		int csize = con.size();
		FOR(i, 0, csize) 
		{
			conf.push_back(p2fea(con[i], general_p));
		}	
		FOR(i, 0, fea.size()) 
		{
			conf[fea[i]].type = fea_p;
		}	
		return conf;
	}


	void protoTile::Trans_contour(vector<Point_f> &c1, Point2f trans_shift)
	{
		for (int i = 0; i < c1.size(); i++)
		{
			c1[i].point+= trans_shift;
		}
	}

	void protoTile::Rotate_contour(vector<Point_f> &c1, Point2f center, double angle)
	{
		vector<Point2f> src = conf_trans(c1);
		vector<Point2f> dst;
		cv::Mat rot_mat = cv::getRotationMatrix2D(center, angle, 1.0);
		/*cv::warpAffine(src, dst, rot_mat, src.size());*/
		cv::transform(src, dst, rot_mat);
		FOR(i, 0, c1.size())
		{
			c1[i].point= dst[i];
		}
	}

	void protoTile::Flip_contour(vector<Point_f> &c1)
	{
		Point2f ccen = center_p(conf_trans(c1));
		int cont_size = c1.size();
		//flip horizontal
		for (int i = 0; i < cont_size; i++)
		{
			c1[i].point.x = 2 * ccen.x - c1[i].point.x;
		}
		for (int i = 0; i < cont_size / 2; i++)
		{
			Point_f mid = c1[i];
			c1[i] = c1[cont_size - 1 - i];
			c1[cont_size - 1 - i] = mid;
		}
	}

	vector<int> protoTile::getCandPoints(vector<Point_f>& contour_f)
	{
		bool show_cand = true;
		vector<int> cand_index;
		int cand_num = 40; //目标候选点个数
		int csize = contour_f.size();
		int margin = csize / cand_num;
		for (int i = 0; i < csize; i += margin) 
		{
			if (contour_f[i].type != fea_p)
			{
				contour_f[i].type = cand_p;
			}
			else continue;
		}
		FOR(i, 0, csize)
		{
			if (contour_f[i].type != general_p)
			{
				cand_index.push_back(i);
			}
		}
		if (show_cand)
		{
			Mat drawing = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
			FOR(i, 0, contour_f.size())
			{
				if (contour_f[i].type == general_p) circle(drawing, contour_f[i].point, 3, Scalar(0, 0, 0), -1);
				else if (contour_f[i].type == fea_p) circle(drawing, contour_f[i].point, 3, Scalar(0, 0, 255), -1);
				else if (contour_f[i].type == cand_p) circle(drawing, contour_f[i].point, 3, Scalar(0, 255, 0), -1);
			}
			imshow("protofirst cand points", drawing);
		}
		return cand_index;
	}
}
