#include "tilingOpt.h"

namespace Tiling_tiles {

	protoTile::protoTile(string filepath, bool show)
	{
		contour = load_point_file(filepath);
		if (contour.empty())
		{
			//cout << filepath << " does not exist! " << endl;
		}
		else
		{
			//cout << "hahahah:" << csize << endl;
			resample(ContourNum, show);
			contour_f = set_flags(contour, feature_points);
		}		
	}

	protoTile::protoTile(vector<Point2f> con)  //根据初始轮廓进行初始化
	{
		contour = con;
		resample(ContourNum);
		contour_f = set_flags(contour, feature_points);
	}

	protoTile::protoTile(vector<Point_f> conf)  //已知点信息的初始化
	{
		contour_f = conf;
		contour = conf_trans(contour_f);
		feature_points = getFeatures(contour_f);
	}


	void protoTile::feature(int  n_min, int n_max, double angle_cos, bool show)
	{
		feature_points = cal_feature(contour, n_min, n_max, angle_cos, show);
	}


	void protoTile::resample(int sam_num, bool show)
	{
		int csize = contour.size();
		feature(1, 0.02*csize, cos(PI * AngleThres / 180), show);
		contour = con_sample(contour, feature_points, sam_num, show);
	}


	void protoTile::scale(double scale)
	{
		Point2f cen = center_p(contour);
		int csize = contour.size();
		for (int i = 0; i < csize; i++)
		{
			contour[i] = scale * (contour[i] - cen) + cen;
		}
		FOR(i, 0, textures.size())
		{
			FOR(j, 0, textures[i].size())
			{
				textures[i][j] = scale * (textures[i][j] - cen) + cen;
				//cout << textures[i][j] << endl;
			}
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


	void protoTile::Trans_proTile(Point2f trans_shift)
	{
		for (int i = 0; i < contour_f.size(); i++)
		{
			contour_f[i].point+= trans_shift;
		}
		for (int i = 0; i < contour.size(); i++)
		{
			contour[i] += trans_shift;
		}
		for (int i = 0; i < textures.size(); i++)
		{
			for (int j = 0; j < textures[i].size(); j++)
			{
				textures[i][j] += trans_shift;
			}
		}		
	}

	void protoTile::Rotate_proTile(Point2f center, double angle)
	{
		cv::Mat rot_mat = cv::getRotationMatrix2D(center, angle, 1.0);
		/*cv::warpAffine(src, dst, rot_mat, src.size());*/
		cv::transform(contour, contour, rot_mat);
		FOR(i, 0, contour_f.size())
		{
			contour_f[i].point= contour[i];
		}
		for (int i = 0; i < textures.size(); i++)
		{
			cv::transform(textures[i], textures[i], rot_mat);
		}
	}

	void protoTile::Flip_proTile()
	{
		Point2f ccen = center_p(contour);
		int cont_size = contour_f.size();
		//flip horizontal
		for (int i = 0; i < cont_size; i++)
		{
			contour_f[i].point.x = 2 * ccen.x - contour_f[i].point.x;
		}
		for (int i = 0; i < cont_size / 2; i++)
		{
			Point_f mid = contour_f[i];
			contour_f[i] = contour_f[cont_size - 1 - i];
			contour_f[cont_size - 1 - i] = mid;
		}

		cont_size = contour.size();
		//flip horizontal
		for (int i = 0; i < cont_size; i++)
		{
			contour[i].x = 2 * ccen.x - contour[i].x;
		}
		for (int i = 0; i < cont_size / 2; i++)
		{
			Point2f mid = contour[i];
			contour[i] = contour[cont_size - 1 - i];
			contour[cont_size - 1 - i] = mid;
		}

		for (int i = 0; i < textures.size(); i++)
		{
			for (int j = 0; j < textures[i].size(); j++)
			{
				textures[i][j].x = 2 * ccen.x - textures[i][j].x;
			}
		}
	}

	vector<int> protoTile::getCandPoints(vector<Point_f>& contour_f)
	{
		bool show_cand = true;
		vector<int> cand_index;
		int cand_num = 30; //目标候选点个数
		int csize = contour_f.size();
		int margin = csize / cand_num;
		for (int i = 0; i < csize; i ++) 
		{
			if (contour_f[i].type != fea_p &&contour_f[max(0, i - 1)].type != fea_p &&contour_f[min(csize - 1, i + 1)].type != fea_p)
			{
				contour_f[i].type = cand_p;
				i = i + margin - 1;
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
		cout << "cand points num: " << cand_index.size() << endl;
		return cand_index;
	}

	vector<int> protoTile::getFeatures(vector<Point_f> &contour_f)
	{
		vector<int> fea_index;
		int csize = contour_f.size();
		FOR(i, 0, csize)
		{
			if (contour_f[i].type == fea_p || contour_f[i].type == fixed_p)
			{
				fea_index.push_back(i);
			}
		}
		return fea_index;
	}

	void protoTile::set_contour(vector<Point2f> c, vector<int> anchor_p, vector<vector<Point2f>> tex)
	{
		contour = c;
		anchor_points = anchor_p;
		textures = tex;
	}

	void protoTile::draw_proTile(Mat &drawing_, Scalar color, Point2f shift, int cut_margin)
	{
		int thickness = 2;
		int lineType = 8;
		int n = contour.size();
		if (cut_margin == 0) 
		{
			for (int t = 0; t < n; t++)
			{
				line(drawing_, contour[t] + shift, contour[(t + 1) % n] + shift, color, thickness, lineType);
			}
		}
		else
		{
			for (int t = 0; t < n; t++)
			{
				Scalar green(0, 255, 100);
				line(drawing_, contour[t] + shift, contour[(t + 1) % n] + shift, green, thickness, lineType);
			}
			for (int anc_in = 0; anc_in < 4; anc_in++)
			{
				int anc_index = anchor_points[anc_in];
				//cout << anc_index << "   " << cut_margin << endl;
				for (int m = 0; m < abs(cut_margin); m++)
				{
					int add_m = (anc_index + cut_margin / abs(cut_margin)*m + n) % n;
					int add_m1 = (anc_index + cut_margin / abs(cut_margin)*(m + 1) + n) % n;
					//cout << add_m<<"   "<< (anc_index + add_m) % n<<"   "<<(anc_index + m) % n << "   " << contour[(anc_index + m) % n] << endl;
					line(drawing_, contour[add_m] + shift, contour[add_m1] + shift, color, thickness, lineType);
				}
			}
		}

		for (int i = 0; i < textures.size(); i++)
		{
			for (int t = 0; t < textures[i].size() -1; t++)
			{
				int lineType = 8;
				line(drawing_, textures[i][t] + shift, textures[i][t + 1] + shift, color, thickness, lineType);
			}
		}
	}

}
