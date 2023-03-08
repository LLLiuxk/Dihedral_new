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
		int all_inner_num = all_inner_conts.size();
		if (all_inner_num == 0)
		{
			std::cout << "no right placement" << endl;
			return;
		}
		FOR(i, 0, 1)
		//FOR(i, 0, all_inner_num)
		{
			match_candidate(i);
			for (int j = 0; j < 1; j++)
			{
				Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
				vector<pair<int, int>> path = cand_paths[i][j];
				vector<Point_f> contour1 = prototile_mid.contour_f;
				vector<Point_f> contour2 = candidate_contours[i][j];
				Point2f sh = Point2f(300, 500) - center_p(conf_trans(contour1));
				draw_pair(draw, conf_trans(contour1), conf_trans(contour2), path, sh);
				
				double con_sc;
				vector<vector<double>> contour2_tar;
				//求contour_对应的point_f和tar值
				//contour2_tar = computeTAR(contour_, con_sc, 0.5);
				double angle = 165;
				double ratio_max = 100;
				double result_score = 0;
				double ratio = 0.5;
				vector<Point_f> contour_2 = morphing_dir(contour1, contour2, path, ratio);
				Point2f sh2 = Point2f(700, 500) - center_p(conf_trans(contour_2));
				draw_contour_points(draw, conf_trans(contour_2), sh2, 3, 2);
				imshow("correspond path", draw);
				vector<Point_f> con_re;
				vector<int> anc1;
				vector<int> anc2;
				FOR(t, 0, contour_2.size())
				{
					if (contour_2[t].type == fixed_p)
						anc1.push_back(t);
				}
				vector<Point2f> frame = { contour_2[anc1[0]].point, contour_2[anc1[1]].point, contour_2[anc1[2]].point, contour_2[anc1[3]].point };
				vector<Point2f> frame_b = base_frame(frame, 3);
				Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
																	  //cout << rot_mat << endl;
				vector<Point2f> contour_dst;
				perspectiveTransform(conf_trans(contour_2), contour_dst, rot_mat);
				
				whole_con_opt(contour_dst, anc1, 0);
				contour_2 = set_flags(contour_dst, contour_2);

				Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
				if (translation_placement(contour_2, con_re, anc1, anc2, draw2))
				{
					cout << "OK!" << endl;
				}

				protoTile c1, c2;
				c1.show_contour(conf_trans(contour_2), anc1);
				c2.show_contour(conf_trans(con_re), anc2);
				RotationVis(c1,c2, AntiClockWise);
				//for (; ratio <= 90; ratio += 5)
				//{
				//	//ratio越高，中间形状保持越多，对应的morphed_A的形状保持越多，eve2越高
				//	double ratio_d = ratio / 100.0;
				//	cout << "ratio_d: " << ratio_d << endl;
				//	//这里final_pettern的点数不一定为200点，但是后续的计算中会重新采样
				//	vector<Point2f> final_pettern = morphing(contour1, contour2, path, ratio_d);
				//	vector<int> return_p;
				//	vector<vector<Point2f>> four_place;
				//	vector<Point2f> morphed_A = extract_contour(final_pettern, mid_inter_morphed, return_p, four_place, all_inner_conts[i].type);
				//	return_p.swap(vector<int>());
				//	four_place.swap(vector<vector<Point2f>>());
				//	vector<Point2f> morphed_A_ori = extract_contour(prototile_mid->contour, mid_inter, return_p, four_place, all_inner_conts[i].type);

				//	double evaluation1 = evalua_deformation(final_pettern, contour_);
				//	double evaluation2 = evalua_deformation(morphed_A, morphed_A_ori);
				//	double score_ave = (1 - ratio_d)*evaluation2 + ratio_d*evaluation1;
				//	if (score_ave > result_score)
				//	{
				//		result_score = score_ave;
				//		ratio_max = ratio_d;
				//	}
				//}

				//vector<Point_f> final_pettern = morphing(contour1, contour2, path, ratio_max);

				//vector<int> return_p;
				//vector<vector<Point2f>> four_place;
				//vector<Point2f> morphed_A = extract_contour(final_pettern, mid_inter_morphed, return_p, four_place, all_inner_conts[i].type);
				////vector<int> return_p2;
				////vector<vector<Point2f>> four_place2;
				////vector<Point2f> morphed_A_ori = extract_contour(prototile_mid->contour, mid_inter, return_p2, four_place2, all_inner_conts[i].type);

				//cout << "ratio_max: " << ratio_max << "  --  result_score:  " << result_score << endl;
				//Mat drawing_pro = Mat(800, 3200, CV_8UC3, Scalar(255, 255, 255));
				//draw_poly(drawing_pro, prototile_mid->contour, Point2f(400, 400));
				//draw_poly(drawing_pro, contour_, Point2f(1200, 400));
				//draw_poly(drawing_pro, final_pettern, Point2f(2000, 400));
				//prototile_tem->loadPoints(final_pettern);
				//draw_poly(drawing_pro, morphed_A, Point2f(2800, 400));
				//string filename = rootname + "\\" + int2string(i) + "_Candidate_" + int2string(j) + "_" + double2string(ratio_max) + "_" + double2string(result_score) + ".png";
				//imwrite(filename, drawing_pro);

				////halftone generate
				//if (!contour_is_simple(morphed_A) || !contour_is_simple(final_pettern))
				//{
				//	cout << "--------------Warning!!!-------------"
				//		<< endl << "The " << j << " candidate patterns are not simple! No results!" << endl << endl;
				//	continue;
				//}
				//int halftone_num = 5000;
				//double l_offset = 2;
				//double scale_t = 0.5;
				//Mat drawing_tesse = Mat(halftone_num, halftone_num, CV_8UC3, Scalar(255, 255, 255));
				//vector<vector<Point2f>> all_tiles = tesse_all(morphed_A, return_p, all_inner_conts[i].type, halftone_num, l_offset, scale_t);
				//int alltilesnum = all_tiles.size();
				//cout << "The num of tiles in final tesse result: " << alltilesnum << endl;
				//for (int m = 0; m < alltilesnum; m++)
				//{
				//	//框架图
				//	//vector<Point2f> tile_dilate = contour_dilate(all_tiles[m], 5);
				//	//vector<Point2f> tile_erode = contour_erode(all_tiles[m], 8);
				//	////draw_poly(drawing_tesse, tile_dilate, center_p(tile_dilate));
				//	////draw_poly(drawing_tesse, tile_erode, center_p(tile_erode), 1);
				//	draw_poly(drawing_tesse, all_tiles[m], center_p(all_tiles[m]));
				//}
				//string filename2 = rootname + "\\" + int2string(i) + "_TessellationResult_" + int2string(j) + ".png";
				//imwrite(filename2, drawing_tesse);
			}
		}


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


	void Tiling_opt::RotationVis(protoTile c1, protoTile c2, int clockorder)
	{
		//       0_______3 0_______3  0________3
		//       |                |  |                 |   |                | 
		//       |      c2      |  |     c1       |   |      c3      | clockwise describes c1's rotating order in opencv
		//       |_______|  |_______|   |_______|
		//       1                2  1               2  1                2
		//
		cout << "c1: "<<c1.contour.size() << c1.contour[0] << " " << c1.contour[1] << " " << c1.contour[2] << " " << c1.contour[3]<<  endl;
		cout << "c2: "<<c2.contour.size() << c2.contour[0] << " " << c2.contour[1] << " " << c2.contour[2] << " " << c2.contour[3] << endl;
		int count = 10;
		while (count > 0)
		{
			Clock_order = clockorder;
			int draw_col = 1000;
			int draw_row = 1000;
			double scale_ratio;
			if (count == 10) scale_ratio = 0.65;
			else scale_ratio = 1.0;
			c1.scale(scale_ratio);
			c2.scale(scale_ratio);
			c1.contour = Trans_contour(c1.contour, Point2f(draw_row / 2, draw_col / 2) - center_p(c1.contour));
			c2.contour = Trans_contour(c2.contour, c1.contour[c1.anchor_points[0]] - c2.contour[c2.anchor_points[3]]);
			//c2.contour = Rotate_contour(c2.contour, center_p(c2.contour), -1);
			if (Clock_order == ClockWise)
			{
				cout << "ClockWise" << endl;
				for (int degree = 0; degree < 46; degree += 2)
				{
					Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
					vector<Point2f> c1_r = Rotate_contour(c1.contour, c1.contour[c1.anchor_points[1]], -degree);
					vector<Point2f> c2_r = Rotate_contour(c2.contour, c2.contour[c2.anchor_points[2]], degree);
					vector<Point2f> c3_r = Trans_contour(c2_r, c1_r[c1.anchor_points[3]] - c2_r[c2.anchor_points[0]]);
					//from top to bottom, from left to right
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[3]] - c1_r[c1.anchor_points[2]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[0]] - c2_r[c2.anchor_points[1]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[3]] - c1_r[c1.anchor_points[2]], 0);
					draw_contour(drawing, c2_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, Point2f(0, 0), 0);
					draw_contour(drawing, c3_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[1]] - c1_r[c1.anchor_points[0]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[2]] - c2_r[c2.anchor_points[3]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[1]] - c1_r[c1.anchor_points[0]], 0);
					imshow("Rotation Visualization", drawing);
					if (degree == 0) waitKey(1000);
					else waitKey(200);
				}
			}
			else if (Clock_order == AntiClockWise)
			{
				cout << "AntiClockWise" << endl;
				for (int degree = 0; degree < 46; degree += 2)
				{
					Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
					vector<Point2f> c1_r = Rotate_contour(c1.contour, c1.contour[c1.anchor_points[0]], degree);
					vector<Point2f> c2_r = Rotate_contour(c2.contour, c2.contour[c2.anchor_points[3]], -degree);
					vector<Point2f> c3_r = Trans_contour(c2_r, c1_r[c1.anchor_points[2]] - c2_r[c2.anchor_points[1]]);
					//from top to bottom, from left to right
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[0]] - c1_r[c1.anchor_points[1]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[3]] - c2_r[c2.anchor_points[2]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[0]] - c1_r[c1.anchor_points[1]], 0);
					draw_contour(drawing, c2_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, Point2f(0, 0), 0);
					draw_contour(drawing, c3_r, Point2f(0, 0), 4);
					draw_contour(drawing, c1_r, c2_r[c2.anchor_points[2]] - c1_r[c1.anchor_points[3]], 0);
					draw_contour(drawing, c2_r, c1_r[c1.anchor_points[1]] - c2_r[c2.anchor_points[0]], 4);
					draw_contour(drawing, c1_r, c3_r[c2.anchor_points[2]] - c1_r[c1.anchor_points[3]], 0);
					//draw_contour(drawing, c2_r, shift2, 3);
					imshow("Rotation Visualization", drawing);
					if (degree == 0) waitKey(1000);
					waitKey(200);
				}
			}
			count--;
		}
		//imwrite("D://Rotation Visualization.png", drawing);
	}


	int Tiling_opt::Tanslation_rule(vector<int> cand_points, vector<Point_f> &contour_s, string rootname)
	{
		int trans = 0;
		int ppindex = cand_points.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(conf_trans(contour_s));
		//std::cout << "contsize: " << cent_cont << endl;
		int times = 0;
		//for (int i = 0; i < ppindex; i++)
		//{
		//	for (int j = i + 1; j < ppindex; j++)
		//	{
		//		for (int m = j + 1; m < ppindex; m++)
		//		{
		//			for (int n = m + 1; n < ppindex; n++)
		//			{					
						int i = 0;
						int j = 8;
						int m = 19;
						int n = 29;
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
						Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (translation_placement( contour_s, inner_contour, indexes, mid_interval, drawing1))
						{
							std::cout << ++trans << " Translation succeed" << endl;
							inPat one_situation(inner_contour, mid_interval, 0);
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(400, 600) - cent_cont;
							FOR(jj, 0, contsize)
							{
								if (contour_s[jj].type == general_p)  circle(drawing1, contour_s[jj].point+ shift2, 2, Scalar(0, 0, 0), -1);
								else if (contour_s[jj].type != general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 255, 0), -1);
							}
							FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);

							string filename = rootname + "/" + to_string(all_inner_conts.size() - 1) + "transPlacingResult.png";
							cv::imwrite(filename, drawing1);
						}
		//			}
		//		}
		//	}
		//}
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
		bool test_coll = true;
		if (test_coll)
		{
			Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
			Point2f cent_cont = center_p(conf_trans(contour_s));
			Point2f shift2 = Point2f(400, 600) - cent_cont;
			draw_contour_points(drawing1, conf_trans(contour_s), shift2);
			FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
			Point2f shift1 = Point2f(1000, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
			FOR(i, 0, 4) draw_poly(drawing1, conf_trans(four_place[i]), shift1);
			string spath = "[" + to_string(indexes[0]) + "," + to_string(indexes[1]) + "," + to_string(indexes[2]) + "," + to_string(indexes[3]) + "]_placement";
			imshow(spath, drawing1);
		}

		if (coll_detec(conf_trans(four_place[0]), conf_trans(four_place[1]), 1) || coll_detec(conf_trans(four_place[0]), conf_trans(four_place[2]), 1)
			|| coll_detec(conf_trans(four_place[0]), conf_trans(four_place[3]), 0) || coll_detec(conf_trans(four_place[1]), conf_trans(four_place[2]), 0))
		{
			//std::cout << "polygon collision" << endl;
			return false;
		}
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
		if (!self_inter)
		{
			FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
			draw_poly(countname, conf_trans(extracted), shift1, 5);
		}
		return !self_inter;
	}


	void Tiling_opt::match_candidate(int inner_index)
	{
		inPat inner= all_inner_conts[inner_index];
		prototile_mid = protoTile(inner.in_contour);
		int cand_num = 10;
		vector<pair<int, bool>> candidate_patterns = compare_TAR(prototile_mid.contour, cand_num);//当前中间图案对应的候选图案
		cout << "This is the " << inner_index << "th inner contour" << endl <<
			"candidate_patterns size: " << candidate_patterns.size() << endl;
		//对中间图案采样200点并重新确认候选点坐标
		//mid_inter = joint_relocate(all_inner_conts[Tiling_index].in_contour, all_inner_conts[Tiling_index].in_interval, num_c);

		double sc_inner = 0;
		vector<vector<double>> inner_tar = compute_TAR(prototile_mid.contour, sc_inner);
		vector<vector<Point_f>> candidate_contour;
		vector<vector<pair<int, int>>> cand_path;

		for (int j = 0; j < cand_num; j++) //只要候选图案里的前cand_num个
		{
			//将所有的结果保存下来
			prototile_second = protoTile(contour_dataset[candidate_patterns[j].first]);
			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//std::cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				prototile_second.Flip_contour(prototile_second.contour_f);
				prototile_second.contour = conf_trans(prototile_second.contour_f);
			}
			else
			{
				//std::cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = match_window_width; //
			double re = tar_mismatch(inner_tar, cand_tar, path, shift, width);

			vector<pair<int, int>> path_min = path;
			/*for (int n = 0; n < path.size(); n++)
			if (prototile_mid->contour_f[path[n].first].type>1)
			path_min.push_back(path[n]);*/

			vector<Point_f> contour_tem_f;
			int c2size = prototile_second.contour_f.size();
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_tem_f.push_back(prototile_second.contour_f[i % c2size]);
			}
			prototile_second.contour_f = contour_tem_f;
			vector<Point2f> contour_cand = conf_trans(prototile_second.contour_f);
			//std::cout << "Mid num: "<<prototile_mid.contour.size() << "   sec num: " << c2size << "  after shifting: " << contour_cand.size() << endl;

			double scale = arcLength(prototile_mid.contour, true) / arcLength(contour_cand, true);
			//std::cout << "scale: " << scale << endl;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
			}
			Point2f cen1 = center_p(prototile_mid.contour);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			double length_min = 1000000;
			int angle_min = 0;
			int times = 0;
			double angl_al = 0;
			Point2f sh_al = Point2f(0, 0);
			Point2f shift_t = Point2f(1000, 1000);
			vector<Point2f> contour_tem;
			//Point2f v0 = contour_inner[0] - cen1;
			//Point2f v1 = contour_cand[0] - cen1;
			//double angle1 = acos(cos_two_vector(v0, v1)) / PI * 180;
			//if (sin_two_vector(v0, v1) < 0) angle1 = -angle1;

			//cout << j << "  path_min_size:  " << path_min.size() << "  cent:" << cen2 << endl;
			while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 100 || angle_min > 10)
			{
				length_min = 1000000;
				angle_min = 0;
				cen2 = center_p(contour_cand);
				//找到距离最近的角度
				for (int angle = 0; angle < 360; angle = angle + 5)
				{
					double leng = 0;
					Mat rot_mat(2, 3, CV_32FC1);
					rot_mat = getRotationMatrix2D(cen2, angle, 1);
					cv::transform(contour_cand, contour_tem, rot_mat);
					for (int m = 0; m < path_min.size(); m++)
					{
						leng += length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]);
					}
					if (leng < length_min)
					{
						angle_min = angle;
						length_min = leng;
					}
				}
				Mat rot_mat1(2, 3, CV_32FC1);
				rot_mat1 = getRotationMatrix2D(cen2, angle_min, 1);
				cv::transform(contour_cand, contour_tem, rot_mat1);
				angl_al += angle_min;
				//contour_cand = contour_mid;
				//移动重心
				shift_t = Point2f(0, 0);
				for (int m = 0; m < path_min.size(); m++)
				{
					shift_t += prototile_mid.contour[path_min[m].first] - contour_tem[path_min[m].second];
				}
				shift_t = Point2f(shift_t.x / path_min.size(), shift_t.y / path_min.size());
				for (int i = 0; i < contour_cand.size(); i++)
				{
					contour_cand[i] = contour_tem[i] + shift_t;
				}
				sh_al += shift_t;
				times++;
			}
			//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << endl;
			prototile_second.contour = contour_cand;
			candidate_contour.push_back(set_flags(contour_cand, prototile_second.contour_f));
			cand_path.push_back(path);
		}
		candidate_contours.push_back(candidate_contour);
		cand_paths.push_back(cand_path);
	}


	vector<pair<int, bool>> Tiling_opt::compare_TAR(vector<Point2f> contour_mid, int chosen_num, int window_width) //chosen_num  选择的最终结果的数目
	{
		//int all_types = 2;
		//std::cout << "inner_c num:" << inner_c.size() << endl;
		double iso_mid = isoperimetric_inequality(contour_mid);
		cout << "Compare with " << all_types << " types candadites." << endl;

		vector<pair<int, bool>> all_final;
		vector<double> all_result;
		double shape_com_mid;
		vector<vector<double>> tar_mid = compute_TAR(contour_mid, shape_com_mid);
		std::cout << "contour_mid: " << contour_mid.size() << "  tar_mid: " << tar_mid.size() << "  each tar: "<< tar_mid[0].size()<<endl;
		//FOR(g, 0, tar_mid[0].size()) cout << tar_mid[0][g] << endl;
		int ttrt = 0;
		for (int index = 0; index < all_types; index++)
		{
			//first filter by isoperimetric_inequality
			double iso_cand = isoperimetric_inequality(conf_trans(contour_dataset[index]));
			double diss = abs(iso_mid - iso_cand) / iso_mid;
			if (diss > SCDisThres)
			{
				ttrt++;
				continue;
			}
			vector<vector<double>> tar_sec = all_con_tars[index];
			vector<vector<double>> tar_sec_f = all_con_tars_flip[index];
			vector<pair<int, int>> path;
			int shift = 0;
			double dis_sc = abs(shape_com_mid - all_shape_complexity[index]);
			//double re = tar_mismatch(tar_mid, tar_sec, mid_index, sec_f, path, shift, window_width);
			double re = tar_mismatch(tar_mid, tar_sec, path, shift, window_width);
			//re = re / (1 + shape_com_mid + all_shape_complexity[index]);  //ori
			re = re / (1 + shape_com_mid + all_shape_complexity[index] - 5 * dis_sc - 2 * diss);
			re = 1 - 0.5*re / path.size();
			//cout << "re: " << re << endl;
			//double re2 = tar_mismatch(tar_mid, tar_sec_f, mid_index, sec_fflip, path, shift, window_width);
			double re2 = tar_mismatch(tar_mid, tar_sec_f, path, shift, window_width);
			//re2 = re2 / (1 + shape_com_mid + all_shape_complexity[index]);
			re2 = re2 / (1 + shape_com_mid + all_shape_complexity[index] - 5 * dis_sc - 2 * diss);
			re2 = 1 - 0.5*re2 / path.size();
			//cout << "re2: " << re2 << endl;
			if (re > re2)
			{
				all_result.push_back(re);
				all_final.push_back(make_pair(index, false));
				//std::cout << re << endl;
			}
			else
			{
				all_result.push_back(re2);
				all_final.push_back(make_pair(index, true));
				//std::cout << re2 << endl;
			}
			//cout << "cand: " << index << "  " << all_result.back() << "  diss:" << diss << endl;
			//system("pause");
		}
		//sort
		double temp;
		pair<int, bool> tempp;
		int all_size = all_result.size();
		for (int i = 0; i < all_size - 1; i++)
			for (int j = 0; j < all_size - 1 - i; j++)
				if (all_result[j] < all_result[j + 1])
				{
					temp = all_result[j];
					all_result[j] = all_result[j + 1];
					all_result[j + 1] = temp;
					tempp = all_final[j];
					all_final[j] = all_final[j + 1];
					all_final[j + 1] = tempp;
				}
		//std::cout << "the fianl order: " << endl;
		vector<pair<int, bool>> all_total_mid;
		int asize = all_final.size();
		chosen_num = min(chosen_num, asize);
		for (int t = 0; t < chosen_num; t++)
		{
			all_total_mid.push_back(all_final[t]);
			std::cout << "order: " << all_final[t].first << "  flip: " << all_final[t].second << " value: " << all_result[t] << " complxeity: " << all_shape_complexity[all_final[t].first] << "     mid sc:" << shape_com_mid << endl;
		}
		cout << "The num of ones byeond ShapeComplexity DisThres: " << ttrt << endl;
		return all_total_mid;
	}

	vector<Point_f> Tiling_opt::morphing(vector<Point_f> contour1, vector<Point_f> contour2, vector<pair<int, int>> path, double ratio)
	{
		//vector<Point2f> final_pettern;
		int cnum1 = contour1.size();
		int cnum2 = contour2.size();
		int psize = path.size();
		//确定以下筛选原则：1.双边选择如果一方有>2的点，则将双方坐标记录下来
		//2.对于多对1的情况：i.如果对面没有>2的，则选择距离最短的点; 2.如果对面有多个，则选择距离最短的; 3如果对面有一个，则选这一个		
		vector<pair<int, int>> final_pair;
		int lock_l = -1;
		int lock_r = -1;  //被锁上的序号不会被再次使用
		for (int j = 0; j < psize - 1; j++)
		{
			//cout << j << "--" << path[j].first << " : " << contour1[path[j].first].type << "    " << path[j].second << " : " << contour2[path[j].second].type << endl;
			if (path[j].first == lock_l || path[j].second == lock_r) continue; //把锁住的情况跳过
			if (contour1[path[j].first].type > 1 || contour2[path[j].second].type > 1)
			{
				int t = j + 1;    //t是下一点
				if (path[t].first != path[j].first && path[t].second != path[j].second)  //两边都没有多对一
				{
					//if (contour1[path[j].first].type > 1) contour2[path[j].second].type = contour1[path[j].first].type;
					//else contour1[path[j].first].type = contour2[path[j].second].type;
					lock_l = path[j].first;
					lock_r = path[j].second;
				}
				else if (path[t].first == path[j].first && path[t].second != path[j].second)  //左边一对多
				{
					double length_min = length_2p(contour1[path[j].first].point, contour2[path[j].second].point);
					int type_max = contour2[path[j].second].type;
					int index = j;
					while (path[t].first == path[j].first)
					{
						if (contour2[path[t].second].type > type_max || (contour2[path[t].second].type == type_max && length_2p(contour1[path[t].first].point, contour2[path[t].second].point) < length_min))
						{
							length_min = length_2p(contour1[path[t].first].point, contour2[path[t].second].point);
							type_max = contour2[path[t].second].type;
							index = t;
						}
						t++;
						if (t == psize) break;
					}
					//lock
					lock_l = path[j].first;
					lock_r = path[index].second;
					j = t - 1;
				}
				else if (path[t].first != path[j].first && path[t].second == path[j].second)  //右边一对多
				{
					double length_min = length_2p(contour1[path[j].first].point, contour2[path[j].second].point);
					int type_max = contour1[path[j].first].type;
					int index = j;
					while (path[t].second == path[j].second)
					{
						if (contour1[path[t].first].type > type_max || (contour1[path[t].first].type == type_max && length_2p(contour1[path[t].first].point, contour2[path[t].second].point) < length_min))
						{
							length_min = length_2p(contour1[path[t].first].point, contour2[path[t].second].point);
							type_max = contour1[path[t].first].type;
							index = t;
						}
						t++;
						if (t == psize) break;
					}
					//lock
					lock_l = path[index].first;
					lock_r = path[j].second;
					j = t - 1;
				}
				if (lock_l != -1 && lock_r != -1) final_pair.push_back(make_pair(lock_l, lock_r));
			}
		}
		int final_pair_size = final_pair.size();

		//设置tolerance，将距离过近的情况合并
		vector<int> de_index;
		for (int mn = 0; mn < final_pair_size; mn++)
		{
			int f_d = 0;
			if (abs((final_pair[mn].first - final_pair[(mn + 1) % final_pair_size].first)) <= tolerance
				|| abs((final_pair[mn].second - final_pair[(mn + 1) % final_pair_size].second)) <= tolerance)
			{
				f_d = 2;
				if (contour1[final_pair[mn].first].type > contour1[final_pair[(mn + 1) % final_pair_size].first].type)
					f_d++;
				else if (contour1[final_pair[mn].first].type < contour1[final_pair[(mn + 1) % final_pair_size].first].type)
					f_d--;
			}
			if (f_d == 3) de_index.push_back((1 + mn++) % final_pair_size);
			if (f_d == 1) de_index.push_back(mn);
			if (f_d == 2 && contour1[final_pair[mn].first].type != 3) de_index.push_back(mn);

		}
		for (int mn = de_index.size() - 1; mn >= 0; mn--)
		{
			pair<int, int> de = delete_vector(final_pair, de_index[mn]);
			cout << de_index[mn] << "  " << de.first << " : " << de.second << endl;
		}
		final_pair_size = final_pair.size();
		cout << "Matched feature point pair : " << final_pair_size << endl;
		//将两个轮廓分段存储
		vector<vector<Point_f>> contour1_seg;
		vector<vector<Point_f>> contour2_seg;

		vector<Point_f> c_seg1;
		vector<Point_f> c_seg2;
		for (int g = 0; g < final_pair_size - 1; g++)
		{
			//cout << final_pair[g].first << "  :  " << final_pair[g].second << endl;
			for (int n = final_pair[g].first; n <= final_pair[g + 1].first; n++)
				c_seg1.push_back(contour1[n]);
			for (int n = final_pair[g].second; n <= final_pair[g + 1].second; n++)
				c_seg2.push_back(contour2[n]);
			contour1_seg.push_back(c_seg1);
			contour2_seg.push_back(c_seg2);
			c_seg1.swap(vector<Point_f>());
			c_seg2.swap(vector<Point_f>());
		}
		for (int n = final_pair[final_pair_size - 1].first; n <= final_pair[0].first + cnum1; n++)
			c_seg1.push_back(contour1[n % cnum1]);
		for (int n = final_pair[final_pair_size - 1].second; n <= final_pair[0].second + cnum2; n++)
			c_seg2.push_back(contour2[n % cnum2]);
		contour1_seg.push_back(c_seg1);
		contour2_seg.push_back(c_seg2);

		//cout << "contour1_seg:  " << contour1_seg.size() << "    contour2_seg:  " << contour2_seg.size() << endl;
		//	<< contour1_seg.back().back().point << "   " << contour1_seg.back().back().type << endl;
		vector<Point_f> contour_fin;
		int ttt = 0;
		contour_fin = morph_segment(contour1_seg[ttt], contour2_seg[ttt], contour1_seg[ttt][0], ratio);
		Mat drawing31 = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
		Point2f shhh = Point2f(600, 600) - contour_fin[0].point;
		line(drawing31, Point2f(600, 600), contour_fin.back().point + shhh, Scalar(50,225,50), 2);
		//for (int nm = 0; nm < contour1_seg[ttt].size(); nm++)
		//	cout << "contour1_seg[g]:" <<contour1_seg[ttt][nm].point << endl;
		//for (int nm = 0; nm < contour2_seg[ttt].size(); nm++)
		//	cout << "contour2_seg[g]:" << contour2_seg[ttt][nm].point << endl;
		//for (int g = 0; g < contour_fin.size(); g++)
		//{
		//	cout << "contour_fin: " << contour_fin[g].point << endl;
		//}
		for (int g = 1; g < contour1_seg.size(); g++)//contour1_seg.size()
		{
			cout << contour1_seg[g].size() << "   :  " << contour2_seg[g].size() << endl << "start before: " << contour_fin.back().point << endl;
			vector<Point_f> contour_tem = morph_segment(contour1_seg[g], contour2_seg[g], contour_fin.back(), ratio);
			line(drawing31, contour_tem[0].point + shhh, contour_tem.back().point + shhh, Scalar(50, 225, 50), 2);
			//for (int nm = 0; nm < contour_tem.size(); nm++)
			//	cout << "contour1_seg[g]:" << contour_tem[nm].point << endl;
			//cout << "end and start : " << contour_tem[0].point << "   :   " << contour_tem.back().point << endl;
			contour_fin.pop_back();
			contour_fin.insert(contour_fin.end(), contour_tem.begin(), contour_tem.end());
		}
		string save = SavePath;
		imwrite(save+"frame.png", drawing31);
		contour_fin.pop_back();

		//final_pettern = p_f2p2f(contour_fin);
		Point2f cen1 = center_p(conf_trans(contour1));
		Point2f cen2 = center_p(conf_trans(contour2));
		//cout << "cen1:  "<<cen1<<"   cen2"<<cen2<<endl;
		//cout << "cnum1:  " << cnum1 << "   cnum2" << cnum1 << "  psize:" << psize<<endl;
		Point2f shift1 = Point2f(400, 400) - cen1;

		Mat drawing3 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		int dd = 0;
		for (int tt = 0; tt < contour1_seg.size(); tt++)//contour1_seg.size()
		{
			//if (tt == 1) continue;
			for (int i = 0; i <contour1_seg[tt].size(); i++)
			{
				//cout << "contour1_seg: "<<contour1_seg[0][i].point + shift1 << endl;
				circle(drawing3, contour1_seg[tt][i].point + shift1, 3, Scalar(255, 0, 0), -1);

				//circle(drawing3, contour_fin[i].point + shift1, 3, Scalar(0, 0, 255), -1);
			}
			for (int i = 0; i <contour2_seg[tt].size(); i++)
			{
				//cout << "contour2_seg: " << contour2_seg[0][i].point + shift1 << endl;
				circle(drawing3, contour2_seg[tt][i].point + shift1, 3, Scalar(0, 255, 0), -1);
			}
		}
		//for (int i = 0; i <contour1_seg[1].size(); i++)
		//{
		//	//cout << "contour1_seg: "<<contour1_seg[0][i].point + shift1 << endl;
		//	circle(drawing3, contour1_seg[1][i].point + shift1, 3, Scalar(155, 120, 0), -1);

		//}
		//for (int i = 0; i <contour2_seg[1].size(); i++)
		//{
		//	//cout << "contour2_seg: " << contour2_seg[0][i].point + shift1 << endl;
		//	circle(drawing3, contour2_seg[1][i].point + shift1, 3, Scalar(0, 155, 120), -1);
		//}
		//MyLine(drawing3, contour1_seg[1][0].point + shift1, contour1_seg[1][0].point + 2*(contour1_seg[1][1].point - contour1_seg[1][0].point) + shift1, "cyan");
		for (int i = 0; i <contour_fin.size(); i++)
		{
			//cout << "contour_fin: " << contour_fin[i].point + shift1 << endl;
			circle(drawing3, contour_fin[i].point + shift1, 3, Scalar(0, 0, 255), -1);
		}
		//MyLine(drawing3, contour2_seg[1][0].point + shift1, contour2_seg[1][0].point + 2 * (contour2_seg[1][1].point - contour2_seg[1][0].point )+ shift1, "cyan");
		imshow("morphing show", drawing3);

		Mat drawing1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i <cnum1; i++)
		{
			if (contour1[i].type == 0)
			{
				circle(drawing1, contour1[i].point + shift1, 3, Scalar(255, 0, 0), -1);
				circle(drawing1, contour2[i].point + shift1, 3, Scalar(255, 125, 0), -1);
			}
			else
			{
				circle(drawing1, contour1[i].point + shift1, 3, Scalar(255, 0, 0), -1);
				circle(drawing1, contour2[i].point + shift1, 3, Scalar(255, 125, 0), -1);
				//circle(drawing1, contour1[i].point + shift1, 5, Scalar(0, 255, 0), -1);
				//circle(drawing1, contour2[i].point + shift1, 5, Scalar(125, 255, 0), -1);
			}
		}
		//MyLine(drawing1, cen1 + shift1, contour1[path[0].first].point + shift1, "red");
		double length = 0;
		for (int j = 0; j < psize; j++)
		{
			if (contour1[path[j].first].type>1)	
				line(drawing1, contour1[path[j].first].point + shift1, contour2[path[j].second].point + shift1, Scalar(100, 100, 100), 2);
			//length += length_two_point2f(contour1[path[j].first].point, contour2[path[j].second].point);
			//cout << path[j].first << " : " << contour1[path[j].first].type << "    " << path[j].second <<" : "<< contour2[path[j].second].type << endl;
			//MyLine(drawing1, contour1[path[j].first].point + shift1, contour2[path[j].second].point + shift1, "gray");
		}

		//cout << "path length: " << length<<endl;
		imshow("path show", drawing1);

		return contour_fin;

	}

	vector<Point_f> Tiling_opt::morph_segment(vector<Point_f> seg1, vector<Point_f> seg2, Point_f start, double ratio) //start 是上一个的尾端，是固定的
	{
		Point_f end = seg1.back();  //如果end.type==3, end 不需要变
		double length_ave = length_2p(start.point, end.point);
		double angle_ave;
		//cout << "1----start: " << start.point << "     --end: " << end.point << endl;
		//确定框架的参数 
		if (end.type != 3) // 确定新的end点 
		{
			double angle1 = cos_2v(Point2f(1, 0), seg1.back().point - seg1[0].point);
			angle1 = (angle1 > 0.999) ? 0.999 : (angle1 < -0.999) ? -0.999 : angle1; //防止出现acos(1)的情况，会返回错误值
			angle1 = acos(angle1);
			if (sin_2v(Point2f(1, 0), seg1.back().point - seg1[0].point) < 0) angle1 = -angle1;
			double angle_12 = cos_2v(seg1.back().point - seg1[0].point, seg2.back().point - seg2[0].point);
			angle_12 = (angle_12 > 0.999) ? 0.999 : (angle_12 < -0.999) ? -0.999 : angle_12;
			angle_12 = acos(angle_12);
			if (sin_2v(seg1.back().point - seg1[0].point, seg2.back().point - seg2[0].point) < 0) angle_12 = -angle_12;
			angle_ave = angle1 + (1 - ratio)*angle_12;
			//angle_ave = angle1 + 0.5*angle_12;
			length_ave = ratio*length_2p(seg1.back().point, seg1[0].point) + (1 - ratio)* length_2p(seg2.back().point, seg2[0].point);
			//cout << "angle1: " << angle1 << " angle_12 " << angle_12 << " length: " << length_ave<<endl;
			Point2f vec_fin = Point2f(cos(angle_ave), sin(angle_ave));
			end.point = start.point + length_ave*vec_fin;
			//end.point = start.point + length_ave*unit_vec(end.point - start.point);
			//end.point = start.point + length_ave*unit_vec(seg1.back().point - seg1[0].point);
		}
		//cout << "2----start: " << start.point << "     --end: " << end.point << endl;
		Point2f mid_des = start.point + 0.5*length_ave * vertical_vec(end.point - start.point);
		vector<Point2f> aff_des;
		aff_des.push_back(start.point);
		aff_des.push_back(end.point);
		aff_des.push_back(mid_des);

		int number_new = 0.5*(seg1.size() + seg2.size());
		//按照number_new进行重采样,这里我们在原样的基础上进行变形
		vector<Point2f> seg_sam1 = sampling_seg(conf_trans(seg1), number_new);
		vector<Point2f> seg_sam2 = sampling_seg(conf_trans(seg2), number_new);
		if (seg_sam1.size() != seg_sam2.size() || seg_sam1.size() != number_new)
		{
			cout << "seg_sam1.size() != seg_sam2.size() != number_new" << endl;
		}
		//x轴上的仿射变化
		Point2f mid_des_x = Point2f(0, 0) + 0.5*length_ave * vertical_vec(Point2f(length_ave, 0));
		vector<Point2f> aff_des_x;
		aff_des_x.push_back(Point2f(0, 0));
		aff_des_x.push_back(Point2f(length_ave, 0));
		aff_des_x.push_back(mid_des_x);
		//计算feature
		/*for (int m = 0; m < 3; m++)
		{
		cout << "aff_des_x: " << aff_des_x[m];
		}*/
		vector<Point2f> aff_seg1;
		aff_seg1.push_back(seg_sam1[0]);
		aff_seg1.push_back(seg_sam1.back());
		aff_seg1.push_back(seg_sam1[0] + 0.5*length_2p(seg_sam1.back(), seg_sam1[0]) * vertical_vec(seg_sam1.back() - seg_sam1[0]));
		Mat M_seg1 = getAffineTransform(aff_seg1, aff_des_x);
		//cout << "M_seg1: " << M_seg1 << endl;
		cv::transform(seg_sam1, seg_sam1, M_seg1);
		vector<Point2f> aff_seg2;
		aff_seg2.push_back(seg_sam2[0]);
		aff_seg2.push_back(seg_sam2.back());
		aff_seg2.push_back(seg_sam2[0] + 0.5*length_2p(seg_sam2.back(), seg_sam2[0]) * vertical_vec(seg_sam2.back() - seg_sam2[0]));
		Mat M_seg2 = getAffineTransform(aff_seg2, aff_des_x);
		cv::transform(seg_sam2, seg_sam2, M_seg2);

		vector<Point2f> contour_morphed;
		for (int j = 0; j < number_new; j++)
		{
			contour_morphed.push_back(ratio*seg_sam1[j] + (1 - ratio)*seg_sam2[j]);
		}
		//计算仿射变换矩阵
		Point2f st_ori = contour_morphed[0];
		Point2f end_ori = contour_morphed.back();
		Point2f mid_ori = st_ori + 0.5*length_2p(st_ori, end_ori)*vertical_vec(end_ori - st_ori);
		vector<Point2f> aff_ori;
		aff_ori.push_back(st_ori);
		aff_ori.push_back(end_ori);
		aff_ori.push_back(mid_ori);
		vector<Point2f> contour_final;
		//将变形得到的feature通过仿射变换移到目的位置
		Mat M1 = getAffineTransform(aff_ori, aff_des);
		cv::transform(contour_morphed, contour_final, M1);
		vector<Point_f> contour_final_;
		FOR(m, 0, contour_final.size())  contour_final_.push_back(p2fea(contour_final[m], 0));
		contour_final_[0].type = seg1[0].type;
		contour_final_.back().type = seg1.back().type;

		return contour_final_;

	}


	vector<Point_f> Tiling_opt::morphing_dir(vector<Point_f> c_mid, vector<Point_f> c_cand, vector<pair<int, int>> path, double ratio)
	{
		Mat draw = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
		draw_contour_points(draw, conf_trans(c_mid), OP, 9);
		draw_contour_points(draw, conf_trans(c_cand), OP, 8);
		cout << "mid_c.size: " << c_mid.size() << "  cand_c.size: " << c_cand.size() << endl;
		int psize = path.size();
		//vector<int> cand_ap = c_mid.anchor_points;// c_cand.anchor_points;
		vector<Point2f> final_mid;
		vector<int> pair_indexs;
		int index = 0; //按照c_mid的顺序
		FOR(i, 0, psize)
		{
			if (path[i].first == index)
			{
				pair_indexs.push_back(path[i].second);
			}
			else
			{
				Point2f new_p;
				/*if (index2 == cand_ap[0] || index2 == cand_ap[1] || index2 == cand_ap[2] || index2 == cand_ap[3])
				{
					new_p = c_cand.contour[index2];
				}
				else
				{
					FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid.contour[pair_indexs[t]] + (1 - ratio) * c_cand.contour[index2]);
					new_p = 1.0 / pair_indexs.size()*new_p;
				}*/
				FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid[index].point + (1 - ratio) * c_cand[pair_indexs[t]].point);
				new_p = 1.0 / pair_indexs.size()*new_p;
				final_mid.push_back(new_p);
				index = path[i].first;
				pair_indexs.swap(vector<int>());
				pair_indexs.push_back(path[i].second);
			}
		}
		if (!pair_indexs.empty())
		{
			Point2f new_p;
			/*if (index2 == cand_ap[0] || index2 == cand_ap[1] || index2 == cand_ap[2] || index2 == cand_ap[3])
			{
				new_p = c_cand.contour[index2];
			}
			else
			{
				FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid.contour[pair_indexs[t]] + (1 - ratio) * c_cand.contour[index2]);
				new_p = 1.0 / pair_indexs.size()*new_p;
			}*/
			FOR(t, 0, pair_indexs.size()) new_p += (ratio *c_mid[index].point + (1 - ratio) * c_cand[pair_indexs[t]].point);
			new_p = 1.0 / pair_indexs.size()*new_p;
			final_mid.push_back(new_p);
		}

		draw_contour_points(draw, final_mid, OP,4);
		imshow("final", draw);
		prototile_fin = protoTile(final_mid);
		prototile_fin.contour_f = set_flags(final_mid, c_mid);

		cout << "final_mid size: " << final_mid.size() << endl;
		return prototile_fin.contour_f;
	}
}