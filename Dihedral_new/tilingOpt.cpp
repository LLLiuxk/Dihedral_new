#include "tilingOpt.h"

namespace Tiling_tiles {

	void Tiling_opt::tiliing_generation(string nameid)
	{
		load_dataset(false);
		string filepath = DefaultPath;
		filepath = filepath + "contour/" + nameid + ".txt";
		prototile_first = protoTile(filepath);
		string savepath = SavePath;
		savepath += nameid;
		const char *na = savepath.c_str();
		if (_access(na, 0) != -1) printf("The  file/dir had been Exisit \n");
		else	_mkdir(na);
		vector<int> cand_points = prototile_first.getCandPoints(prototile_first.contour_f);
		int trans = Tanslation_rule(cand_points, prototile_first.contour_f, savepath);
		//int flips = Flipping_rule(p_p_index, cont_orig, rootname);


	}

	void Tiling_opt::load_dataset(bool input_images)
	{
		if (!contour_dataset.empty() || !all_con_tars.empty())
		{
			std::cout << "the Contours OR TARs have been exist" << endl;
			contour_dataset.swap(vector<vector<Point_f>>());
			all_con_tars.swap(vector<vector<vector<double>>>());
			all_con_tars_flip.swap(vector<vector<vector<double>>>());
			all_shape_complexity.swap(vector<double>());
		}
		if (input_images)
		{
			//vector<int> t = { 26 };
			//int tin = 0;
			FOR(i, 0, AllTypes)
			{			
				//if (i != t[tin]) continue;
				//else tin++;
				string imagepath = DefaultPath;
				imagepath = imagepath + "dataset/" + to_string(i) + ".png";
				ima2contour(imagepath, input_images);
			}
		}
		FOR(i, 0, AllTypes)
		{
			string filepath = DefaultPath;
			filepath = filepath + "contour/" + to_string(i) + ".txt";
			cout << "Image " << to_string(i) << endl;
			protoTile read_tile(filepath);
			if (read_tile.contour.empty()) continue;
			contour_dataset.push_back(read_tile.contour_f);
			double shape_com;
			vector<vector<double>> tar_all = compute_TAR(read_tile.contour, shape_com);
			vector<vector<double>> tar_all_flip = compute_TAR(Flip_contour(read_tile.contour), shape_com);
			all_con_tars.push_back(tar_all);
			all_con_tars_flip.push_back(tar_all_flip);
			all_shape_complexity.push_back(shape_com);
		}
		all_types = contour_dataset.size();
		std::cout << "All types： " << all_types << "  load over!" << endl;
		std::cout << "All TARs of contours have been computed" << endl;
	}

	int Tiling_opt::Tanslation_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname)
	{
		int trans = 0;
		int ppindex = cand_points.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(conf_trans(contour_s));
		//std::cout << "contsize: " << cent_cont << endl;
		int times = 0;
		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				for (int m = j + 1; m < ppindex; m++)
				{
					for (int n = m + 1; n < ppindex; n++)
					{
						times++;
						//std::cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
						//if (abs(part_points_index[n] - part_points_index[m]) < margin) continue;
						vector<Point_f> inner_contour;
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						vector<int> indexes; //1.translati on：以下所有该类摆放都以1-3,2-4为轴摆放
						indexes.push_back(cand_points[i]);
						indexes.push_back(cand_points[j]);
						indexes.push_back(cand_points[m]);
						indexes.push_back(cand_points[n]);
						Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (translation_placement( contour_s, inner_contour, indexes, mid_interval, drawing1))
						{
							std::cout << ++trans << " Translation succeed" << endl;
							allnum_inner_c++;
							inPat one_situation(inner_contour, mid_interval, 0);
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(400, 400) - cent_cont;
							FOR(jj, 0, contsize)
							{
								if (contour_s[i].type == general_p)  circle(drawing1, contour_s[i].point+ shift2, 2, Scalar(0, 0, 0), -1);
								else if (contour_s[i].type != general_p)  circle(drawing1, contour_s[i].point + shift2, 2, Scalar(0, 255, 0), -1);
							}
							FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);

							string filename = rootname + "/" + to_string(allnum_inner_c - 1) + "transPlacingResult.png";
							cv::imwrite(filename, drawing1);
						}
					}
				}
			}
		}
		std::cout << "Trans times: " << times << endl;
		return trans;
	}

	int Tiling_opt::Flipping_rule(vector<int> part_points_index, vector<Point_f> &contour_s, string rootname)
	{

	}

	bool Tiling_opt::translation_placement(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname)
	{
		int csize = contour_s.size();
		Point2f line1 = contour_s[indexes[2]].point - contour_s[indexes[0]].point;
		Point2f line2 = contour_s[indexes[3]].point - contour_s[indexes[1]].point;
		vector<Point2f> shifting;
		shifting.push_back(line1);
		shifting.push_back(line2);
		shifting.push_back(line1 + line2);
		// 提取平移围成的轮廓
		vector<vector<Point_f>> four_place;
		four_place.push_back(contour_s);
		for (int i = 0; i < 3; i++)
		{
			vector<Point_f> one_loca;
			for (int j = 0; j < csize; j++)
			{
				one_loca.push_back(Point_f(contour_s[j].point + shifting[i], contour_s[j].type));
			}
			four_place.push_back(one_loca);
		}
		//std::cout << "compute four_place" << endl;
		//translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
		//              ①     ③            共有①②③④四个图案，
		//               /\      /\			   检测①中的顶点3与②中顶点1，
		//             /1 \   /1 \           以及①中顶点4与③中顶点2的角度之和
		//            /      \/      \         
		//          2\   4/\2    /4
		//              \  /    \  /
		//             3\/      \/3
		//               /\      /\     
		//             /1 \   /1 \ 
		//            /     \/       \          为保证在一定范围内可靠，
		//          2\   4/\2  4/          比较时计算∠Pt-iPtPt+i的角度，
		//              \  /   \  /           i从1取到3
		//              3\/    \/3
		//               ②    ④
		int total_num = 0;
		ex_indexes.swap(vector<int>());
		extracted.swap(vector<Point_f>());
		ex_indexes.push_back(0);
		for (int t = indexes[3]; t > indexes[2]; t--)                   //translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
		{
			extracted.push_back(four_place[0][t]);
			total_num++;
			//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[0] + csize; t > indexes[3]; t--)
		{
			extracted.push_back(four_place[1][t % csize]);
			total_num++;
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[1]; t > indexes[0]; t--)
		{
			extracted.push_back(four_place[3][t]);
			total_num++;
		}
		ex_indexes.push_back(total_num);
		for (int t = indexes[2]; t > indexes[1]; t--)
		{
			extracted.push_back(four_place[2][t]);
		}
		//提取后检测是否自交
		int index_1, index_2;
		bool self_inter = self_intersect(conf_trans(extracted), index_1, index_2);
		Point2f shift1 = Point2f(1000, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
		FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
		draw_poly(countname, conf_trans(extracted), shift1, 5);
		string spath = "[" + to_string(indexes[0]) + "," + to_string(indexes[1]) + "," + to_string(indexes[2]) + "," + to_string(indexes[3]) + "]_placement";
		if (self_inter)
		{
			//cout << "intersection!" << endl;	
			//imshow(spath, countname);
		}
	
		return !self_inter;
	}

}