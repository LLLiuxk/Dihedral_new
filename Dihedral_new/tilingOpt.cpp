#include "tilingOpt.h"
vector<string> frame_type={ "parallelogram", "rhombus", "rectangle", "square" };
bool pers_trans = 0;
bool coll_opt = 0;
bool deve_opt = 0;

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
		std::cout << "---------------------------------------" << endl << "Matching Candadites" << endl << "---------------------------------------" << endl;
		FOR(i, 0, all_inner_num)
		{
			match_candidate(i, savepath);
			//for (int j = 0; j < 5; j++)
			//{
			//	vector<pair<int, int>> path = cand_paths[j];
			//	vector<Point_f> contour1 = prototile_mid.contour_f;
			//	vector<Point_f> contour2 = candidate_contours[j];

			//	prototile_second = protoTile(contour2);
			//	vector<pair<int, int>> path_fea = cand_fea_paths[j];

			//	//show the path
			//	draw_path(prototile_mid.contour, prototile_second.contour, prototile_mid.feature_points, prototile_second.feature_points, path_fea);
			//	
			//	double con_sc;
			//	vector<vector<double>> contour2_tar;
			//	//求contour_对应的point_f和tar值
			//	//contour2_tar = computeTAR(contour_, con_sc, 0.5);
			//	double angle = 165;
			//	double ratio_max = 100;
			//	double result_score = 0;
			//	double ratio = 0.5;
			//	double min_score = 1000;
			//	vector<Point_f> contour_2;
			//	vector<Point_f> con_re;
			//	vector<int> anc_mid;
			//	vector<int> anc_re;

			//}
		}


	}



	void Tiling_opt::tiliing_gen_specify(string nameid, vector<int> anc_points)
	{
		load_dataset(false);
		string filepath = DefaultPath;
		filepath = filepath + "contour/" + nameid + ".txt";
		prototile_first = protoTile(filepath);
		string savepath = SaveSpecPath;
		savepath = savepath +nameid+"/";
		const char *na = savepath.c_str();
		if (_access(na, 0) != -1) printf("The  file/dir had been Exisit \n");
		else	_mkdir(na);
		vector<int> cand_points = prototile_first.getCandPoints(prototile_first.contour_f);
		int trans = Tanslation_rule_spec(cand_points, prototile_first.contour_f, savepath, anc_points);
		//int flips = Flipping_rule(p_p_index, cont_orig, rootname);
		int all_inner_num = all_inner_conts.size();
		if (all_inner_num == 0)
		{
			std::cout << "no right placement" << endl;
			return;
		}
		std::cout << "---------------------------------------" << endl << "Matching Candadites" << endl << "---------------------------------------" << endl;
		FOR(i, 0, 1)
			//FOR(i, 0, all_inner_num)
		{
			match_candidate_(i);
			//feature_match();
			for (int j = 0; j < 1; j++)
			{
				vector<pair<int, int>> path = cand_paths[j];
				vector<Point_f> contour1 = prototile_mid.contour_f;
				vector<Point_f> contour2 = candidate_contours[j];
				prototile_second = protoTile(contour2);
				vector<pair<int, int>> path_fea = cand_fea_paths[j];

				//show the path
				draw_path(prototile_mid.contour, prototile_second.contour, prototile_mid.feature_points, prototile_second.feature_points, path_fea);

				double con_sc;
				vector<vector<double>> contour2_tar;
				//求contour_对应的point_f和tar值
				//contour2_tar = computeTAR(contour_, con_sc, 0.5);
				double angle = 165;
				double ratio_max = 100;
				double result_score = 0;
				double ratio = 0.5;
				double min_score = 1000;
				vector<Point_f> contour_2;
				vector<Point_f> con_re;
				vector<int> anc_mid;
				vector<int> anc_re;
				Mat draw_min;

				for (double r = 0.1; r <= 0.91; r += 0.05) 
				{
					//ratio越高，中间形状保持越多，对应的morphed_A的形状保持越多，eve2越高
					cout << "ratio: " << r << "   ";
					vector<Point_f> c_2 = morphing(contour1, contour2, path_fea, r);

					//calculate two deformed shapes
					vector<Point_f> c_re;
					vector<int> a_mid;
					vector<int> a_re;
					int csize = c_2.size();
					FOR(t, 0, csize)
					{
						if (c_2[t].type == fixed_p)
							a_mid.push_back(t);
					}
					Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
					if (!translation_spec(c_2, c_re, a_mid, a_re, draw2))
					{
						cout << "Warning! Intersection Occur!" << endl;
						//continue;
					}
					/*Mat overlap_map = Mat(800, 1800, CV_8UC3, Scalar(255, 255, 255));
					draw_overlap(overlap_map, conf_trans(c_2), conf_trans(contour1), Point2f(0,0));
					draw_overlap(overlap_map, conf_trans(c_2), conf_trans(contour2), Point2f(600, 0));
					draw_overlap(overlap_map, conf_trans(contour1), conf_trans(contour2), Point2f(1200, 0));
					imwrite("D:/vs2015project/Dihedral_new/maps/" + to_string(r) + "over.png", overlap_map);
					cout<< evaluation_area(conf_trans(contour1), conf_trans(contour2))<<"    "<< evaluation_area(conf_trans(c_2), conf_trans(contour1))<<"    "<< evaluation_area(conf_trans(c_2), conf_trans(contour2))<<endl;*/
					double score1 = deform_evalue(conf_trans(c_2), conf_trans(contour1));
					double score2 = deform_evalue(conf_trans(c_2), conf_trans(contour2));
					//double score = sqrt(pow(score1, 2) + pow(score2, 2));
					double score = max(score1, score2);
					cout << "score: "<< score << endl;
					if (min_score > score)
					{
						min_score = score;
						ratio = r;
						contour_2 = c_2;
						con_re = c_re;
						anc_mid = a_mid;
						anc_re = a_re;
						draw_min = draw2;
					}
				}
				imshow( "[" + to_string(anc_mid[0]) + "," + to_string(anc_mid[1]) + "," + to_string(anc_mid[2]) + "," + to_string(anc_mid[3]) + "]_placement", draw_min);

				cout << endl << "min ratio: " << ratio << "   min score: " << min_score << "   c1&c2 size:  " << contour1.size() << "  " << contour2.size() << "  contour_2.size: " << contour_2.size() << endl << endl;
				Mat pair_match = Mat(600, 1000, CV_8UC3, Scalar(255, 255, 255));
				Point2f sh = Point2f(300, 300) - center_p(conf_trans(contour1));
				draw_pair(pair_match, conf_trans(contour1), conf_trans(contour2), path, sh);
				Point2f sh2 = Point2f(700, 300) - center_p(conf_trans(contour_2));
				draw_contour_points(pair_match, conf_trans(contour_2), sh2, 3, 2);
				imshow("correspond path", pair_match);
				Mat two_c = Mat(1000, 1600, CV_8UC3, Scalar(255, 255, 255));
				Point2f ima_sh = Point2f(400, 500) - center_p(conf_trans(contour_2));
				draw_contour_points(two_c, conf_trans(contour_2), ima_sh, 5, 2);
				draw_contour_points(two_c, conf_trans(con_re), ima_sh, 6, 2);

				//optimization for two contours
				std::cout << "---------------------------------------" << endl << "Optimize The Dihedral Tessellation" << endl << "---------------------------------------" << endl;
				Clock_order = AntiClockWise;  //指示初始input形状旋转的方向
				int Clock_order2 = abs(Clock_order - 1);
				con_re = contour_opt(con_re, anc_re, 1, 0, 1, 0, 0, Clock_order);
				contour_2 = contour_opt(contour_2, anc_mid, 1, 1, 1, 0, 0, Clock_order2);
				//cout << "size: " << contour_2.size() << "  " << con_re.size() << endl;
				//将两个轮廓对齐
				Point2f shift = contour_2[anc_mid[3]].point - con_re[anc_re[0]].point;
				FOR(ii, 0, contour_2.size()) con_re[ii].point += shift;

				FOR(gg, 0, 4)  cout << contour_2[anc_mid[gg]].point << "    " << con_re[anc_re[gg]].point << endl;
				ima_sh = Point2f(1200, 500) - center_p(conf_trans(contour_2));
				draw_contour_points(two_c, conf_trans(contour_2), ima_sh, 5, 2);
				FOR(ii, 0, 4) draw_contour_points(two_c, conf_trans(con_re), ima_sh + contour_2[anc_mid[ii]].point - con_re[anc_re[(ii + 1) % 4]].point, 6, 2);
				imshow("2 contours:", two_c);

				double merge_ratio = 0.5;
				double min_mers = 1000;
				for (double mr = 0.2; mr < 0.81; mr += 0.05)
				{
					vector<Point_f> cont2_new = contour_2;
					vector<Point_f> conre_new = con_re;
					vector<int> anc1 = anc_mid;
					vector<int> anc2 = anc_re;
					//merge_contours(cont2_new, conre_new, anc1, anc2, mr);
					merge_contours(conre_new, cont2_new, anc2, anc1,  mr);
					double merge_s1 = deform_evalue(conf_trans(cont2_new), conf_trans(contour_2));
					double merge_s2 = deform_evalue(conf_trans(conre_new), conf_trans(con_re));
					//double merge_score = sqrt(pow(merge_s1, 2) + pow(merge_s2, 2));
					double merge_score = max(merge_s1, merge_s2);
					if (merge_score < min_mers) 
					{
						min_mers = merge_score;
						merge_ratio = mr;
					}
					cout << "mr: "<<mr<<"   merge_score: " << merge_score << endl;
				}
				cout << endl << "min_ratio: " << merge_ratio << "   " << min_mers << endl << endl;

				string file_name ="Two_contour.svg";
				std::ofstream file(file_name);
				file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";
				FOR(p_index, 0, con_re.size())
				{
					file << "<line x1=\"" << con_re[p_index].point.x << "\" y1=\"" << con_re[p_index].point.y << "\" x2=\"" << con_re[(p_index + 1)% con_re.size()].point.x << "\" y2=\"" << con_re[(p_index + 1) % con_re.size()].point.y << "\" stroke=\"rgb(120,120,120)\" />\n";
					if(con_re[p_index].type==fixed_p)
						file << "<circle cx=\"" << con_re[p_index].point.x << "\" cy=\"" << con_re[p_index].point.y << "\" r=\"2\" fill=\"rgb(240, 176, 30)\" />\n";
				}
				FOR(p_index, 0, contour_2.size())
				{
					file << "<line x1=\"" << contour_2[p_index].point.x - 300 << "\" y1=\"" << contour_2[p_index].point.y << "\" x2=\"" << contour_2[(p_index + 1) % contour_2.size()].point.x - 300 << "\" y2=\"" << contour_2[(p_index + 1) % contour_2.size()].point.y << "\" stroke=\"rgb(120,120,120)\" />\n";
					if (contour_2[p_index].type == fixed_p)
						file << "<circle cx=\"" << contour_2[p_index].point.x - 300 << "\" cy=\"" << contour_2[p_index].point.y << "\" r=\"2\" fill=\"rgb(240, 176, 30)\" />\n";
				}
				file << "</svg>";
				file.close();
				//merge_contours(contour_2, con_re, anc_mid, anc_re, merge_ratio);
				merge_contours(con_re, contour_2,  anc_re, anc_mid, merge_ratio);
				//FOR(gg, 0, 4)  cout << contour_2[anc_mid[gg]].point << "    " << con_re[anc_re[gg]].point << endl;
				
				file_name = "Merged_contour.svg";
				std::ofstream file1(file_name);
				file1 << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";
				FOR(p_index, 0, con_re.size())
				{
					file1 << "<line x1=\"" << con_re[p_index].point.x << "\" y1=\"" << con_re[p_index].point.y << "\" x2=\"" << con_re[(p_index + 1) % con_re.size()].point.x << "\" y2=\"" << con_re[(p_index + 1) % con_re.size()].point.y << "\" stroke=\"rgb(120,120,120)\" />\n";
					if (con_re[p_index].type == fixed_p)
						file1 << "<circle cx=\"" << con_re[p_index].point.x << "\" cy=\"" << con_re[p_index].point.y << "\" r=\"2\" fill=\"rgb(240, 176, 30)\" />\n";
				}
				FOR(p_index, 0, contour_2.size())
				{
					file1 << "<line x1=\"" << contour_2[p_index].point.x - 300 << "\" y1=\"" << contour_2[p_index].point.y << "\" x2=\"" << contour_2[(p_index + 1) % contour_2.size()].point.x - 300 << "\" y2=\"" << contour_2[(p_index + 1) % contour_2.size()].point.y << "\" stroke=\"rgb(120,120,120)\" />\n";
					if (contour_2[p_index].type == fixed_p)
						file1 << "<circle cx=\"" << contour_2[p_index].point.x - 300 << "\" cy=\"" << contour_2[p_index].point.y << "\" r=\"2\" fill=\"rgb(240, 176, 30)\" />\n";
				}
				file1 << "</svg>";
				file1.close();

				write_twoCon(savepath + "twoCon.txt", anc_re, conf_trans(con_re), anc_mid, conf_trans(contour_2));
				string command = "DDET.exe  " + to_string(Clock_order) + "  "+ savepath;
				system(command.c_str());
				vector<Point2f> con1 = conf_trans(con_re);
				vector<Point2f> con2 = conf_trans(contour_2);
				con1 = load_point_file(savepath + "c1.txt");
				con2 = load_point_file(savepath + "c2.txt");
				if (con1.empty() || con2.empty())
				{
					cout << "No contours! Exit!" << endl;
					exit(0);
				}
				vector<vector<Point2f>> tex1 = load_texture(savepath + "texture1.txt");
				vector<vector<Point2f>> tex2 = load_texture(savepath + "texture2.txt");
				contour2obj(savepath + "c1.obj", con1, Point2f(200, 200) - center_p(con1), 15);
				contour2obj(savepath + "c2.obj", con2, Point2f(200, 200) - center_p(con1), 15);

				protoTile c1, c2;
				c1.set_contour(con1, anc_re, tex1);
				c2.set_contour(con2, anc_mid, tex2);
				RotationVis(c1, c2, Clock_order, savepath);

			}
		}
	}


	vector<Point_f> Tiling_opt::contour_opt(vector<Point_f> cont, vector<int>& anc_p, int type, int times, bool pers_trans , bool coll_opt, bool deve_opt, int cworder) //type: 0=contours bbx; 1: square bbx
	{
		vector<Point_f> con_re;
		vector<Point2f> contour_dst = conf_trans(cont);
		int csize = contour_dst.size();

		vector<vector<int>> handle_area;
		vector<vector<Point2f>> handle_points;
		FOR(index_num, 0, anc_p.size())
		{
			int area_width = handle_area_width;
			vector<int> area_one;
			vector<Point2f> handle_one;
			FOR(width_n, -area_width, area_width + 1) 
			//for (int width_n = -area_width; width_n < area_width + 1; width_n += 2)
			{
				int index_ = (anc_p[index_num] + width_n + csize) % csize;
				area_one.push_back(index_);
				//handle_points.push_back(contour_2[index_].point + shift_new);
				handle_one.push_back(cont[index_].point);
			}
			handle_area.push_back(area_one);
			handle_points.push_back(handle_one);
		}
		MatrixXd V;
		MatrixXi F;
		//vector<Point2f> tt = triangulate_Contours_bbx(cont_re, anc_re, V, F);
		//vector<Point2f> tt = triangulate_Contours_bbx(contour_dst, anc_mid, V, F);
		vector<Point2f> tt;
		if (type == 0)
			tt = triangulateContour(contour_dst, V, F);
		else if(type==1)
			tt = triangulate_bbx(contour_dst, V, F);
		else if (type == 2)
			tt= triangulate_Contours_bbx(contour_dst, anc_p, V, F);
	/*	else if (type == 3)
		{
			tt = triangulate_2Contours(vector<Point2f>& cont1, vector<Point2f>& cont2, MatrixXd& V, MatrixXi& F);
		}*/
		Mat image = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
		Point stf=Point2f(600, 600) - center_p(contour_dst);
		for (size_t i = 0; i < F.rows(); i++)
		{
			Point pt1(V.row(F(i, 0)).x(), V.row(F(i, 0)).y());
			Point pt2(V.row(F(i, 1)).x(), V.row(F(i, 1)).y());
			Point pt3(V.row(F(i, 2)).x(), V.row(F(i, 2)).y());
			//cout << F(i, 0) << "   " << F(i, 1) << "   " << F(i, 2) << endl;
			//cout << pt1 << "    " << pt2 << "   " << pt3 << endl;
			line(image, pt1 + stf, pt2 + stf, Scalar(0, 180, 0), 1, LINE_AA);
			line(image, pt2 + stf, pt3 + stf, Scalar(0, 180, 0), 1, LINE_AA);
			line(image, pt3 + stf, pt1 + stf, Scalar(0, 180, 0), 1, LINE_AA);
		}
		//imwrite("Triangulation.png", image);
		imshow("Triangulation", image);
		double degree_after_opt = 0;
		if (pers_trans)
		{
			cout << "Triangulation and ARAP deforming......" << endl;
			vector<Point2f> frame = { contour_dst[anc_p[0]], contour_dst[anc_p[1]], contour_dst[anc_p[2]], contour_dst[anc_p[3]] };
			//vector<Point2f> frame = { cont_re[anc_re[0]], cont_re[anc_re[1]], cont_re[anc_re[2]], cont_re[anc_re[3]] };
			if (type == 2)
			{
				Point2f sh1 = frame[3] - frame[0];
				FOR(n, 0, 4) frame[n] += sh1;
			}
			int min_type = 0;
			double min_l = 10000;
			vector<Point2f> frame_b = base_frame(frame, min_type);
			FOR(f_type, 0, 4)
			{
				vector<Point2f> frame_bb = base_frame(frame, f_type);
				vector<pair<int, int>> path_min = { make_pair(0,0),make_pair(1,1),make_pair(2,2),make_pair(3,3) };
				double alige_e = conotour_align(frame, frame_bb, path_min);
				if (alige_e < min_l)
				{
					min_l = alige_e;
					min_type = f_type;
					frame_b = frame_bb;
				}
			}
			cout << "frame_type: " << frame_type[min_type] << "   " << min_l << endl;
			vector<int> anc_mid_;
			FOR(g, 0, 4)
			{
				//cout << anc_p[g]<<"   "<<contour_dst[anc_p[g]] << "   " << tt[anc_p[g]] << "   " << frame[g] << endl;
				anc_mid_.push_back(point_locate(tt, frame[g]));
				cout << "frame[g]: " << frame[g] <<"   "<< anc_mid_.back()<< endl;
			}
			cout << frame[1] - frame[0] << "   " << frame[2] - frame[1] << endl;
			// before 
			Mat draw_ = Mat(1000, 1000, CV_8UC3, Scalar(255, 255, 255));
			Point2f sh1 = Point2f(500, 500) - center_p(contour_dst);
			draw_contour(draw_, contour_dst, sh1, 5);
			FOR(m, 0, 4) circle(draw_, frame[m] + sh1, 3, Scalar(125, 0, 0));

			string parap = ParaPath;
			string obj_path = parap + "mid_result.obj";
			string para_path = parap + "deform_para.txt";
			string deformed_c = parap + "deformed_c.txt";
			write_obj(obj_path, V, F);
			write_para(para_path, anc_mid_, frame_b);
			//write_para(para_path, anc_mid, frame_b);
			//write_para(para_path, handle_area, handle_points);
			string command = "D:/vs2015project/Dihedral_new/Dihedral_new/ARAP_Deform.exe  " + to_string(times) + "  " + obj_path + "  " + para_path + " " + deformed_c+"  "+to_string(csize);
			//string command = "D:/vs2015project/ARAP_Deform/x64/Debug/ARAP_Deform.exe  " + obj_path + "  " + para_path + " " + deformed_c;
			cout << command << endl;
			system(command.c_str());
			contour_dst = load_point_file(deformed_c);
			vector<Point2f> con_tem;
			FOR(m, 0, csize) con_tem.push_back(contour_dst[m]);
			contour_dst = con_tem;
			//recover the boundary
			string file_name = to_string(times) + "_After_align.svg";
			std::ofstream file(file_name);
			file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";
			for (int i = 0; i < csize; i++)
			{
				file << "<line x1=\"" << contour_dst[i].x << "\" y1=\"" << contour_dst[i].y << "\" x2=\"" << contour_dst[(i + 1) % csize].x << "\" y2=\"" << contour_dst[(i + 1) % csize].y << "\" stroke=\"rgb(120,120,120)\" />\n";
				file << "<circle cx=\"" << contour_dst[i].x << "\" cy=\"" << contour_dst[i].y << "\" r=\"2\" fill=\"rgb(100, 100, 100)\" />\n";
			}
			FOR(index_h, 0, handle_area.size())
			{
				vector<int> new_indexs = handle_area[index_h];
				vector<Point2f> new_handp;
				//cout << "index_h: " << index_h << "   " << handle_area.size() << endl;
				FOR(i,0, new_indexs.size())
				{
					new_handp.push_back(contour_dst[new_indexs[i]]);
					//cout << "i: " << i << "   " << contour_dst[new_indes[i]] << endl;
				}
				
				bound_recover(handle_points[index_h], new_handp);

				FOR(i, 0, new_indexs.size())
					contour_dst[new_indexs[i]] = new_handp[i];
				file << "<line x1=\"" << contour_dst[(new_indexs[0] - 1 + csize)%csize].x << "\" y1=\"" << contour_dst[(new_indexs[0] - 1 + csize) % csize].y << "\" x2=\"" << contour_dst[new_indexs[0]].x << "\" y2=\"" << contour_dst[new_indexs[0]].y << "\" stroke=\"rgb(149, 202, 135)\" />\n";
				FOR(i, 0, new_indexs.size())
				{
					file << "<line x1=\"" << contour_dst[new_indexs[i]].x << "\" y1=\"" << contour_dst[new_indexs[i]].y << "\" x2=\"" << contour_dst[(new_indexs[i] + 1) % csize].x << "\" y2=\"" << contour_dst[(new_indexs[i] + 1) % csize].y << "\" stroke=\"rgb(149, 202, 135)\" />\n";
					file << "<circle cx=\"" << contour_dst[new_indexs[i]].x << "\" cy=\"" << contour_dst[new_indexs[i]].y << "\" r=\"2\" fill=\"rgb(250, 65, 65)\" />\n";
				}
				/*FOR(i, 0, new_indexs.size())
					file << "<circle cx=\"" << contour_dst[new_indexs[i]].x << "\" cy=\"" << contour_dst[new_indexs[i]].y << "\" r=\"2\" fill=\"rgb(250, 65, 65)\" />\n";*/
			}
			file << "</svg>";
			file.close();
			//resample contour_dst
			
			Mat draw_22 = Mat(800, 1200, CV_8UC3, Scalar(255, 255, 255));
			vector<Point2f> resam_ = { contour_dst[anc_p[0]],contour_dst[anc_p[1]], contour_dst[anc_p[2]], contour_dst[anc_p[3]] };
			Point2f sht = Point2f(300, 300) - center_p(contour_dst);
			draw_contour_points(draw_22, contour_dst,sht);
			draw_contour_points(draw_22, resam_, sht, 5,3);
			circle(draw_22, contour_dst[0]+sht, 3, Scalar(0, 255, 0), 1);
			protoTile mid(contour_dst);
			contour_dst = mid.contour;
			vector<int> anc_re = relocate(resam_, contour_dst);
			int move_sht = 0 - anc_re[0];
			if (move_sht != 0)
			{
				vector<Point2f> new_c;
				FOR(ii, 0, csize) new_c.push_back(contour_dst[(anc_re[0] + ii) % csize]);
				contour_dst = new_c;
				FOR(ii, 0, 4) anc_re[ii] = (anc_re[ii] + move_sht + csize) % csize;
			}
			sht = Point2f(600, 0) + sht;
			draw_contour_points(draw_22, contour_dst, sht);
			circle(draw_22, contour_dst[0] + sht, 3, Scalar(0, 255, 0), 1);
			FOR(tt,0,4) circle(draw_22, contour_dst[anc_re[tt]] + sht, 3, Scalar(0, 0, 255), 3);
			//imwrite("D:/ressample.png", draw_22);
			FOR(ii, 0, 4)
			{
				cout << anc_p[ii]<<"   "<<resam_[ii] << "   " << anc_re[ii]<<"   "<< contour_dst[anc_p[ii]]<<"   "<<contour_dst[anc_re[ii]] << endl;
				contour_dst[anc_re[ii]] = resam_[ii];
			}
			anc_p = anc_re;
			FOR(ii, 0, csize) con_re.push_back(Point_f(contour_dst[ii], general_p));
			FOR(jj, 0, anc_re.size()) con_re[anc_re[jj]].type = fixed_p;

			cout << "csize " << csize << "  "<<contour_dst.size() <<"  "<<con_re.size()<< endl;
			draw_contour(draw_, contour_dst, sh1, 8);
			FOR(m, 0, 4) circle(draw_, frame_b[m] + sh1, 3, Scalar(0, 255, 0));
			imshow(to_string(times)+" After align handle points:", draw_);
			//con_re = set_flags(contour_dst, cont);
		}
		if (coll_opt)
		{
			cout << "contour fine tuning......" << endl;
			//contour_de_crossing(contour_dst);
			string file_name = to_string(times) + "_Fine_tuning.svg";
			std::ofstream file(file_name);
			file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";
			for (int i = 0; i < csize; i++)
			{
				file << "<line x1=\"" << contour_dst[i].x << "\" y1=\"" << contour_dst[i].y << "\" x2=\"" << contour_dst[(i + 1) % csize].x << "\" y2=\"" << contour_dst[(i + 1) % csize].y << "\" stroke=\"rgb(120,120,120)\" />\n";
				file << "<circle cx=\"" << contour_dst[i].x << "\" cy=\"" << contour_dst[i].y << "\" r=\"2\" fill=\"rgb(100, 100, 100)\" />\n";
			}

			contour_fine_tuning(contour_dst);
			//con_re = set_flags(contour_dst, cont);
			con_re = set_flags(contour_dst, con_re);

			for (int i = 0; i < csize; i++)
			{
				file << "<line x1=\"" << contour_dst[i].x << "\" y1=\"" << contour_dst[i].y << "\" x2=\"" << contour_dst[(i + 1) % csize].x << "\" y2=\"" << contour_dst[(i + 1) % csize].y << "\" stroke=\"rgb(100,100,200)\" />\n";
				file << "<circle cx=\"" << contour_dst[i].x << "\" cy=\"" << contour_dst[i].y << "\" r=\"2\" fill=\"rgb(100, 100, 100)\" />\n";
			}

			cout << "csize " << csize << "  " << contour_dst.size() << "  " << con_re.size() << endl;
		}
		if (deve_opt)
		{
			cout << "contour deployability optimizing......" << endl;
			//degree_after_opt = whole_con_opt(resam_dst, anc_re, 0);
			degree_after_opt = whole_con_opt(contour_dst, anc_p, cworder);
			//con_re = set_flags(contour_dst, cont);
			con_re = set_flags(contour_dst, con_re);
			Mat dep_opt  = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
			draw_contour_points(dep_opt, contour_dst, Point2f(400, 400) - center_p(contour_dst), 5, 2);
			FOR(m, 0, 4) circle(dep_opt, contour_dst[anc_p[m]] + Point2f(400, 400) - center_p(contour_dst), 3, Scalar(0, 255, 0));
			imshow(to_string(times) + " After contour deployability optimization:", dep_opt);
			//degree_after_opt = whole_con_opt(contour_dst, anc_mid, 0);
			cout << endl<<"After developable optimation, the collision degree: " << degree_after_opt << endl;
			cout << "csize " << csize << "  " << contour_dst.size() << "  " << con_re.size() << endl;
		}
		return con_re;
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
			FOR(i, 171, 172)
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
			cout << "Image " << to_string(i) << "   ";
			protoTile read_tile(filepath, false);
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

	void Tiling_opt::load_para(std::string filename)
	{
		ifstream fin(filename);
		if (fin.fail() == true)
		{
			cout << "Cannot open load_pare file." << endl;
		}
		string param;
		while (fin >> param)
		{
			if (param == "image_id")
				fin >> image_id;
		}
		cout << "Paramaters Loaded Successfully." << endl << endl;
	}


	void Tiling_opt::RotationVis(protoTile c1, protoTile c2, int clockorder, string save_path)
	{
		//       0_______3 0_______3  0________3
		//       |                |  |                 |   |                | 
		//       |      c2      |  |     c1       |   |      c3      | clockwise describes c1's rotating order in opencv
		//       |_______|  |_______|   |_______|
		//       1                2  1               2  1                2
		//
		cout << "c1: "<<c1.contour.size() << "   "<<c1.contour[c1.anchor_points[0]] << "  " << c1.contour[c1.anchor_points[1]]
			<< " " << c1.contour[c1.anchor_points[2]] << "   " << c1.contour[c1.anchor_points[3]]<<  endl;
		cout << "c2: "<<c2.contour.size() << "   " << c2.contour[c2.anchor_points[0]] << " " << c2.contour[c2.anchor_points[1]] 
			<< " " << c2.contour[c2.anchor_points[2]] << " " << c2.contour[c2.anchor_points[3]] << endl;
		int rot_deg = 46;
		vector<Mat> images;
		int draw_col = 1000;
		int draw_row = 1000;
		double scale_ratio = 0.65;
		c1.scale(scale_ratio);
		c2.scale(scale_ratio);
		//对齐到图案中央
		c1.Trans_proTile(Point2f(draw_row / 2, draw_col / 2) - center_p(c1.contour));
		c2.Trans_proTile(c1.contour[c1.anchor_points[0]] - c2.contour[c2.anchor_points[3]]);

		int add_degree = 1;
		int cut_mar = CutMargin;
		vector<double> Angles;
		vector<double> Prs;
		vector<Point2f> c1_total_ori;

		if (clockorder == ClockWise)
		{
			cout << "ClockWise" << endl;
			for (int degree = 0; degree < rot_deg; degree += add_degree)
			{
				if (degree != 0) cut_mar = 0;
				Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
				protoTile c3 = c2;
				c3.Trans_proTile(c1.contour[c1.anchor_points[3]] - c2.contour[c2.anchor_points[0]]);			
				//from top to bottom, from left to right
				vector<Point2f> shifts = { 
					c2.contour[c2.anchor_points[3]] - c1.contour[c1.anchor_points[2]] ,
					c1.contour[c1.anchor_points[0]] - c2.contour[c2.anchor_points[1]],
					c3.contour[c2.anchor_points[3]] - c1.contour[c1.anchor_points[2]],
					OP,OP,OP,
					c2.contour[c2.anchor_points[1]] - c1.contour[c1.anchor_points[0]],
					c1.contour[c1.anchor_points[2]] - c2.contour[c2.anchor_points[3]],
					c3.contour[c2.anchor_points[1]] - c1.contour[c1.anchor_points[0]]
				};
				c1.draw_proTile(drawing, colorbar[0].second, shifts[0], -cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[1], cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[2], -cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[3], cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[4], -cut_mar);
				c3.draw_proTile(drawing, colorbar[4].second, shifts[5], cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[6], -cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[7], cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[8], -cut_mar);

				if (degree % 2 == 0)
				{
					vector<Point2f> c1_total = mergedVector({ c1.contour,c1.contour,c1.contour,c1.contour }, { shifts[0],shifts[2],shifts[6],shifts[8] });
					if (degree == 0) c1_total_ori = c1_total;
					else
					{
						double pr = cal_Poisson_ratio(c1_total_ori, c1_total);
						Angles.push_back(degree);
						Prs.push_back(pr);
						//cout << "DEGREE: " << degree << "    " << pr << endl;
					}
				}

				c1.Rotate_proTile(c1.contour[c1.anchor_points[1]], -add_degree);
				c2.Rotate_proTile(c2.contour[c2.anchor_points[2]], add_degree);
				
				imshow("Rotation Visualization", drawing);
				images.push_back(drawing);
				if (degree == 0) waitKey(800);
				else waitKey(100);
			}
		}
		else if (clockorder == AntiClockWise)
		{
			cout << "AntiClockWise" << endl;
			for (int degree = 0; degree < rot_deg; degree += add_degree)
			{
				if (degree != 0) cut_mar = 0;
				Mat drawing = Mat(draw_row, draw_col, CV_8UC3, Scalar(255, 255, 255));
				protoTile c3 = c2;
				c3.Trans_proTile(c1.contour[c1.anchor_points[2]] - c2.contour[c2.anchor_points[1]]);
				//from top to bottom, from left to right
				vector<Point2f> shifts = {
					c2.contour[c2.anchor_points[0]] - c1.contour[c1.anchor_points[1]],
					c1.contour[c1.anchor_points[3]] - c2.contour[c2.anchor_points[2]],
					c3.contour[c2.anchor_points[0]] - c1.contour[c1.anchor_points[1]],
					OP,OP,OP,
					c2.contour[c2.anchor_points[2]] - c1.contour[c1.anchor_points[3]],
					c1.contour[c1.anchor_points[1]] - c2.contour[c2.anchor_points[0]],
					c3.contour[c2.anchor_points[2]] - c1.contour[c1.anchor_points[3]]
				};
				c1.draw_proTile(drawing, colorbar[0].second, shifts[0], cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[1], -cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[2], cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[3], -cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[4], cut_mar);
				c3.draw_proTile(drawing, colorbar[4].second, shifts[5], -cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[6], cut_mar);
				c2.draw_proTile(drawing, colorbar[4].second, shifts[7], -cut_mar);
				c1.draw_proTile(drawing, colorbar[0].second, shifts[8], cut_mar);

				/*if(degree)*/
				if (degree % 2 == 0)
				{
					vector<Point2f> c1_total = mergedVector({ c1.contour,c1.contour,c1.contour,c1.contour }, { shifts[0],shifts[2],shifts[6],shifts[8] });
					if (degree == 0) c1_total_ori = c1_total;
					else
					{
						double pr = cal_Poisson_ratio(c1_total_ori, c1_total);
						Angles.push_back(degree);
						Prs.push_back(pr);
						//cout << "DEGREE: " << degree << "    " << pr << endl;
					}
				}

				c1.Rotate_proTile(c1.contour[c1.anchor_points[0]], add_degree);
				c2.Rotate_proTile(c2.contour[c2.anchor_points[3]], -add_degree);
				
				imshow("Rotation Visualization", drawing);
				images.push_back(drawing);
				if (degree == 0) waitKey(800);
				waitKey(100);
			}
		}

		string file_name = "result.svg";
		std::ofstream file1(file_name);
		file1 << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";
		file1 << "<polygon points=\"";

		for (const auto& point : c1.contour) {
			file1 << point.x << "," << point.y << " ";
		}
		file1 << "\" stroke=\"black\"  fill=\"rgb(199,144,172)\" />\n";
		file1 << "<polygon points=\"";

		for (const auto& point : c2.contour) {
			file1 << point.x << "," << point.y << " ";
		}
		file1 << "\" stroke=\"black\"  fill=\"rgb(202,225,195)\" />\n";
		file1 << "</svg>";
		file1.close();


		drawLineGraph(Angles, Prs, save_path + "NPR_angle.png");
		write_avi(images, save_path + "output.avi", 5);
		imwrite(save_path + "Rotation Visualization.png", images[0]);
		imshow("Rotation Visualization2222", images[0]);

		Mat drawing_all = Mat(5000, 5000, CV_8UC3, Scalar(255, 255, 255));
		//draw_repeat(drawing_all, vector<Point2f> c1, vector<int> anc1, vector<Point2f> c2, vector<int> anc2);


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
			progress_bar(i, ppindex);
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
						Mat drawing1 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (translation_placement( contour_s, inner_contour, indexes, mid_interval, drawing1))
						{
							//std::cout << ++trans << " Translation succeed" << endl;
							inPat one_situation(inner_contour, mid_interval, 0);
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(300, 600) - cent_cont;
							FOR(jj, 0, contsize)
							{
								if (contour_s[jj].type == general_p)  circle(drawing1, contour_s[jj].point+ shift2, 2, Scalar(0, 0, 0), -1);
								else if (contour_s[jj].type != general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 255, 0), -1);
							}
							FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
							circle(drawing1, contour_s[0].point + shift2, 4, Scalar(255, 0, 55), -1);

							string filename = rootname + "/" + to_string(all_inner_conts.size() - 1) + "transPlacingResult.png";
							cv::imwrite(filename, drawing1);
						}
					}
				}
			}
		}
		std::cout << "Trans times: " << times << endl;
		return trans;
	}

	int Tiling_opt::Tanslation_rule_spec(vector<int> cand_points, vector<Point_f> &contour_s, string rootname, vector<int> anc_points)
	{
		int trans = 0;
		int ppindex = cand_points.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(conf_trans(contour_s));
		//std::cout << "contsize: " << cent_cont << endl;
		int times = 0;
		int i = anc_points[0];
		int j = anc_points[1];
		int m = anc_points[2];
		int n = anc_points[3];
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

		if (translation_placement(contour_s, inner_contour, indexes, mid_interval, drawing1))
		{
			//std::cout << ++trans << " Translation succeed" << endl;
			inPat one_situation(inner_contour, mid_interval, 0);
			all_inner_conts.push_back(one_situation);

			Point2f shift2 = Point2f(300, 600) - cent_cont;
			FOR(jj, 0, contsize)
			{
				if (contour_s[jj].type == general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 0, 0), -1);
				else if (contour_s[jj].type != general_p)  circle(drawing1, contour_s[jj].point + shift2, 2, Scalar(0, 255, 0), -1);
			}
			FOR(jj, 0, 4) circle(drawing1, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
			circle(drawing1, contour_s[0].point + shift2, 4, Scalar(255, 0, 55), -1);

			string filename = rootname + to_string(all_inner_conts.size() - 1) + "transPlacingResult.png";
			cv::imwrite(filename, drawing1);
			cv::imshow("Initial tiling placement: ", drawing1);
		}
		std::cout << "Trans times: " << times << endl;
		return trans;
	}


	int Tiling_opt::Flipping_rule(vector<int> part_points_index, vector<Point_f> &contour_s, string rootname)
	{
		return 0;
	}

	bool Tiling_opt::translation_placement(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname)
	{
		bool test_coll = false;
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
		Point2f shift1 = Point2f(1100, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
		if (!self_inter)
		{
			FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
			draw_poly(countname, conf_trans(extracted), shift1, 5);
		}
		return !self_inter;
	}

	bool Tiling_opt::translation_spec(vector<Point_f> &contour_s, vector<Point_f> &extracted, vector<int> indexes, vector<int> &ex_indexes, Mat &countname)
	{
		int csize = contour_s.size();
		//cout << "contour_s size: " << csize << endl;
		Point2f line1 = contour_s[indexes[2]].point - contour_s[indexes[0]].point;
		Point2f line2 = contour_s[indexes[3]].point - contour_s[indexes[1]].point;
		vector<Point2f> shifting;
		shifting.push_back(line1);
		shifting.push_back(line1 + line2);
		shifting.push_back(line2);

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
		Point2f cent_cont = center_p(conf_trans(contour_s));
		Point2f shift2 = Point2f(400, 600) - cent_cont;
		draw_contour_points(countname, conf_trans(contour_s), shift2);
		FOR(jj, 0, 4) circle(countname, contour_s[indexes[jj]].point + shift2, 3, Scalar(0, 0, 255), -1);
		Point2f shift1 = Point2f(1000, 600) - (center_p(conf_trans(contour_s)) + 0.5*(line1 + line2));
		FOR(i, 0, 4) draw_poly(countname, conf_trans(four_place[i]), shift1);
		string spath = "[" + to_string(indexes[0]) + "," + to_string(indexes[1]) + "," + to_string(indexes[2]) + "," + to_string(indexes[3]) + "]_placement";

		int total_num = 0;
		ex_indexes.swap(vector<int>());
		extracted.swap(vector<Point_f>());
		//ex_indexes.push_back(0);
		for (int i = 0; i < four_place.size(); i++)
		{
			ex_indexes.push_back(total_num);
			int s_i = (i + 3) % 4;
			int e_i = (i + 2) % 4;
			int t = indexes[s_i];
			if (s_i<e_i)  t += csize;
			for (; t > indexes[e_i]; t--)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
			{
				extracted.push_back(four_place[i][t % csize]);
				total_num++;
				//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
			}
		}
		//for (int t = indexes[3]; t > indexes[2]; t--)                   //translation的提取规律为1的4-3, 2的1-4, 3的2-1, 4的3-2
		//{
		//	extracted.push_back(four_place[0][t]);
		//	total_num++;
		//	//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[0] + csize; t > indexes[3]; t--)
		//{
		//	extracted.push_back(four_place[1][t % csize]);
		//	total_num++;
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[1]; t > indexes[0]; t--)
		//{
		//	extracted.push_back(four_place[3][t]);
		//	total_num++;
		//}
		//ex_indexes.push_back(total_num);
		//for (int t = indexes[2]; t > indexes[1]; t--)
		//{
		//	extracted.push_back(four_place[2][t]);
		//}
		draw_poly(countname, conf_trans(extracted), shift1, 5);
		imshow(spath, countname);
		//提取后检测是否自交
		int index_1, index_2;
		bool self_inter = self_intersect(conf_trans(extracted), index_1, index_2);
		return !self_inter;
	}

	void Tiling_opt::match_candidate(int inner_index, string rootname)
	{
		inPat inner = all_inner_conts[inner_index];
		prototile_mid = protoTile(inner.in_contour);
		int cand_num = 10;
		int cal_cand_num = 5;
		vector<pair<int, bool>> candidate_patterns = compare_TAR(prototile_mid.contour_f, cand_num);//当前中间图案对应的候选图案
		cout << "This is the " << inner_index << "th inner contour" << endl <<
			"candidate_patterns size: " << candidate_patterns.size() << endl;

		double sc_inner = 0;
		vector<vector<double>> inner_tar = compute_TAR(prototile_mid.contour, sc_inner);
		//vector<vector<Point_f>> candidate_contour;
		//vector<vector<pair<int, int>>> cand_path;
		candidate_contours.swap(vector<vector<Point_f>>());
		cand_paths.swap(vector<vector<pair<int, int>>>());
		cand_fea_paths.swap(vector<vector<pair<int, int>>>());
		path_shift.swap(vector<int>());

		Mat draw_cands = Mat(600, 3000, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(draw_cands, conf_trans(inner.in_contour), Point2f(250, 300) - center_p(conf_trans(inner.in_contour)), 5);
		cal_cand_num = min(cal_cand_num, (int)candidate_patterns.size());
		for (int j = 0; j < cal_cand_num; j++) //只要候选图案里的前cal_cand_num个
		{
			//将所有的结果保存下来
			prototile_second = protoTile(contour_dataset[candidate_patterns[j].first]);

			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//std::cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				prototile_second.Flip_proTile();
				prototile_second.feature_points = prototile_second.getFeatures(prototile_second.contour_f);
			}
			else
			{
				//std::cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
			}
			int cmsize = prototile_mid.contour_f.size();
			vector<int> mid_index(cmsize, 0), sec_f(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n] = prototile_mid.contour_f[n].type;
				sec_f[n] = prototile_second.contour_f[n].type;
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = match_window_width; //
			double re = tar_mismatch_fea(inner_tar, cand_tar, mid_index, sec_f, path, shift, width);
			vector<pair<int, int>> path_min = path;
			vector<Point_f> contour_tem_f;
			int c2size = prototile_second.contour_f.size();
			//移动点顺序
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_tem_f.push_back(prototile_second.contour_f[i % c2size]);
			}
			prototile_second.contour_f = contour_tem_f;

			vector<Point2f> contour_cand = conf_trans(prototile_second.contour_f);
			//std::cout << "Mid num: "<<prototile_mid.contour.size() << "   sec num: " << c2size << "  after shifting: " << contour_cand.size() << endl;

			//移动到重心重叠的地方
			Point2f cen1 = center_p(prototile_mid.contour);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			//two types of calculating the scale: i.based on the arclength  ii.based in uniform scale
			double scale = arcLength(prototile_mid.contour, true) / arcLength(contour_cand, true);
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
			}

			double length_min = 1000000;
			int angle_min = 0;
			int times = 0;
			double angl_al = 0;
			Point2f sh_al = Point2f(0, 0);
			Point2f shift_t = Point2f(1000, 1000);
			vector<Point2f> contour_tem;
			//cout << j << "  path_min_size:  " << path_min.size() << "  cent:" << cen2 << endl;
			while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 10 || angle_min > 10)// 
			{
				length_min = 1000000;
				angle_min = 0;
				cen2 = center_p(contour_cand);
				//找到距离最近的角度
				for (int angle = 0; angle < 360; angle = angle + 2)
				{
					double leng = 0;
					Mat rot_mat(2, 3, CV_32FC1);
					rot_mat = getRotationMatrix2D(cen2, angle, 1);
					cv::transform(contour_cand, contour_tem, rot_mat);
					for (int m = 0; m < path_min.size(); m++)
					{
						leng += pow(length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]), 2);
						/*leng += length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]);*/
					}
					leng = sqrt(leng);
					if (leng < length_min)
					{
						angle_min = angle;
						length_min = leng;
					}
				}
				//cout << "angle_min: " << angle_min << "   length_min:  " << length_min << endl;
				Mat rot_mat1(2, 3, CV_32FC1);
				rot_mat1 = getRotationMatrix2D(cen2, angle_min, 1);
				cv::transform(contour_cand, contour_tem, rot_mat1);
				angl_al += angle_min;
				//contour_cand = contour_tem;
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
				//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << "   shift_t: " << shift_t << endl;
			}
			//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << endl;
			prototile_second.contour = contour_cand;
			prototile_second.contour_f = set_flags(contour_cand, prototile_second.contour_f);
			candidate_contours.push_back(prototile_second.contour_f);
			cand_paths.push_back(path);

			cout << "The " << j << "th cands' path order:" << endl;
			//filter out inappropriate matches
			filter_path(path, prototile_second.contour_f, prototile_mid.contour_f, inner_tar, cand_tar, shift);

			Point2f shift_cand = Point2f(500 * j + 750, 300)-center_p(contour_cand);
			draw_poly(draw_cands, contour_cand, shift_cand);
		}
		string filename = rootname + "/cands/";
		_mkdir(filename.c_str());
		filename = filename + to_string(inner_index) + "_candidates.png";
		cv::imwrite(filename, draw_cands);
	}

	void Tiling_opt::match_candidate_(int inner_index)
	{
		inPat inner = all_inner_conts[inner_index];
		prototile_mid = protoTile(inner.in_contour);
		int cand_num = 10;
		vector<pair<int, bool>> candidate_patterns = compare_TAR(prototile_mid.contour_f, cand_num);//当前中间图案对应的候选图案
		cout << "This is the " << inner_index << "th inner contour" << endl <<
			"candidate_patterns size: " << candidate_patterns.size() << endl;

		double sc_inner = 0;
		vector<vector<double>> inner_tar = compute_TAR(prototile_mid.contour, sc_inner);
		//vector<vector<Point_f>> candidate_contour;
		//vector<vector<pair<int, int>>> cand_path;
		candidate_contours.swap(vector<vector<Point_f>>());
		cand_paths.swap(vector<vector<pair<int, int>>>());
		cand_fea_paths.swap(vector<vector<pair<int, int>>>());
		path_shift.swap(vector<int>());

		for (int j = 0; j < 1; j++) //只要候选图案里的前cand_num个
		{
			//将所有的结果保存下来
			prototile_second = protoTile(contour_dataset[candidate_patterns[j].first]);
			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//std::cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				prototile_second.Flip_proTile();
				prototile_second.feature_points = prototile_second.getFeatures(prototile_second.contour_f);
			}
			else
			{
				//std::cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
			}
			int cmsize = prototile_mid.contour_f.size();
			vector<int> mid_index(cmsize, 0), sec_f(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n] = prototile_mid.contour_f[n].type;
				sec_f[n] = prototile_second.contour_f[n].type;
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = match_window_width; //
			double re = tar_mismatch_fea(inner_tar, cand_tar, mid_index, sec_f, path, shift, width); 
			vector<pair<int, int>> path_min = path;
			vector<Point_f> contour_tem_f;
			int c2size = prototile_second.contour_f.size();
			//移动点顺序
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_tem_f.push_back(prototile_second.contour_f[i % c2size]);
			}
			prototile_second.contour_f = contour_tem_f;

			vector<Point2f> contour_cand = conf_trans(prototile_second.contour_f);
			//std::cout << "Mid num: "<<prototile_mid.contour.size() << "   sec num: " << c2size << "  after shifting: " << contour_cand.size() << endl;

			//移动到重心重叠的地方
			Point2f cen1 = center_p(prototile_mid.contour);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			//two types of calculating the scale: i.based on the arclength  ii.based in uniform scale
			double scale = arcLength(prototile_mid.contour, true) / arcLength(contour_cand, true);
			//std::cout << "scale: " << scale << endl;
			/*double s_1 = 0, s_2 = 0;
			FOR(mn, 0, contour_cand.size())
			{
				s_1 += pow(length_2p(prototile_mid.contour[mn], cen1), 2);
				s_2 += pow(length_2p(contour_cand[mn], cen1), 2);
			}
			double scale2 = sqrt(s_1 / s_2);
			cout << "scale2: " << scale2 << endl;*/
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
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
			while (times < 3 || length_2p(shift_t, Point2f(0, 0)) > 10 || angle_min > 10)// 
			{
				length_min = 1000000;
				angle_min = 0;
				cen2 = center_p(contour_cand);
				//找到距离最近的角度
				for (int angle = 0; angle < 360; angle = angle + 2)
				{
					double leng = 0;
					Mat rot_mat(2, 3, CV_32FC1);
					rot_mat = getRotationMatrix2D(cen2, angle, 1);
					cv::transform(contour_cand, contour_tem, rot_mat);
					for (int m = 0; m < path_min.size(); m++)
					{
						leng += pow(length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]),2);
						/*leng += length_2p(prototile_mid.contour[path_min[m].first], contour_tem[path_min[m].second]);*/
					}
					leng = sqrt(leng);
					if (leng < length_min)
					{
						angle_min = angle;
						length_min = leng;
					}
				}
				//cout << "angle_min: " << angle_min << "   length_min:  " << length_min << endl;
				Mat rot_mat1(2, 3, CV_32FC1);
				rot_mat1 = getRotationMatrix2D(cen2, angle_min, 1);
				cv::transform(contour_cand, contour_tem, rot_mat1);
				angl_al += angle_min;
				//contour_cand = contour_tem;
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
				cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << "   shift_t: "<< shift_t<<endl;
			}
			//cout << "angle_al: " << angl_al << "   --  " << length_min << "    sh_al: " << sh_al << endl;
			prototile_second.contour = contour_cand;
			prototile_second.contour_f = set_flags(contour_cand, prototile_second.contour_f);
			candidate_contours.push_back(prototile_second.contour_f);
			cand_paths.push_back(path);

			//filter out inappropriate matches
			filter_path(path, prototile_second.contour_f, prototile_mid.contour_f, inner_tar, cand_tar, shift);
			
		}
	}


	void Tiling_opt::filter_path(vector<pair<int, int>> path, vector<Point_f> con2, vector<Point_f> con_mid, vector<vector<double>> inner_tar, vector<vector<double>> cand_tar, int shift)
	{
		/*筛选策略：
		i. 遍历path，如果1中是fea_p，将last_index到范围内的特征点加到waiting队列，如果1中是fixed_p，将范围内的所有特征点以及path对应点加到waiting中
		ii. 从waiting中选出匹配cost最小的一个，与path_fea中对比，如果1是fixed，直接把path_fea中的交叉或者重叠路径删除，如果是fea，必不会交叉，但是如果有重叠，需计算cost*/
		vector<pair<int, int>> path_fea;
		int path_size = path.size();
		int sec_size = con2.size();
		//cout << path_size << "   " << sec_size << endl;
		int mar = path_margin;
		int last_index = -mar;
		int last_m = -1;
		FOR(m, 0, path_size)
		{
			//if (m > 0) last_first = path[m - 1].first;
			//cout << "m: " << m << "   " << path[m].first << "  " << path[m].second << endl
			//	<< "   type:  "<<prototile_mid.contour_f[path[m].first].type << "  " << prototile_second.contour_f[path[m].second].type << endl;
			pair<int, int> one_pair = path[m];
			int sec_index = one_pair.second;
			if (con_mid[one_pair.first].type == fixed_p || con_mid[one_pair.first].type == fea_p)  //must have a match
			{
				vector<int> waiting_merge;
				//cout << " last_m:  " << last_m << "  last_index  " << last_index << endl;
				if (abs(last_m - last_index) > mar) last_index = last_index - (last_m + sec_size);
				//cout << "last_index  " << last_index << endl;
				FOR(n, 1 - mar, mar)
				{
					int tem_sec = (sec_index + n + sec_size) % sec_size;
					if (con_mid[one_pair.first].type == fea_p)
					{
						if (tem_sec >= last_index && con2[tem_sec].type == fea_p)
							waiting_merge.push_back(tem_sec);
					}	
					if (con_mid[one_pair.first].type == fixed_p)
					{
						if (tem_sec == sec_index || con2[tem_sec].type == fea_p)
							waiting_merge.push_back(tem_sec);
					}
				}
				//if (m == 85) cout << waiting_merge.size()<<"    "<<waiting_merge[0]  << endl;
				//the waiting_merge of fixed_p must not be empty
				if (!waiting_merge.empty())
				{
					vector<Point2f> con_mid_c = conf_trans(con_mid);
					vector<Point2f> con2_c = conf_trans(con2);
					//FOR(mm, 0, waiting_merge.size()) cout << one_pair.first << "  waiting_merge[mm]: " << waiting_merge[mm] << endl;
					int pcsize = con_mid_c.size();
					int min_index = 0;
					if (waiting_merge.size() > 1)
					{
						double min_dis_tar = 10000;
						FOR(t, 0, waiting_merge.size())
						{
							//calculate the pair with less tar distance
							double normal_r = 0.5 / (2 * PI);
							int cand_index = (waiting_merge[t] + shift) % cand_tar.size();
							double ang1 = angle_2v(con_mid_c[(one_pair.first - 1 + pcsize) % pcsize] - con_mid_c[one_pair.first], con_mid_c[(one_pair.first + 1) % pcsize] - con_mid_c[one_pair.first]);
							double ang2 = angle_2v(con2_c[(waiting_merge[t] - 1 + pcsize) % pcsize] - con2_c[waiting_merge[t]], con2_c[(waiting_merge[t] + 1) % pcsize] - con2_c[waiting_merge[t]]);
							//double dis_tar = tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]) + 0.02*abs(dis_cos1 - dis_cos2);
							double diss_angle = abs(ang1 - ang2);
							if (ang1*ang2 < 0) diss_angle = 2 * PI - diss_angle;  //<0, 异号
							double dis_tar = tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]) + normal_r*diss_angle;
							//cout << "distar: " << dis_tar << "   " << tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]) << "  ang1: " << ang1 / PI * 180 << "    " << ang2 / PI * 180 << "   " << normal_r*diss_angle << endl;
							if (dis_tar < min_dis_tar)
							{
								min_dis_tar = dis_tar;
								min_index = t;
							}
						}
					}
					//cout << "fea m: " << m << "   min_index: " << min_index << "    " << waiting_merge[min_index] << endl;
					if (con_mid[one_pair.first].type == fixed_p)
					{
						int wait_chosen = waiting_merge[min_index];
						if (!path_fea.empty() &&path_fea.back().second >= waiting_merge[min_index] && abs(path_fea.back().second - waiting_merge[min_index])<2*path_margin)
						{
							if (con_mid[path_fea.back().first].type == fixed_p)  wait_chosen = wait_chosen + (one_pair.first - path_fea.back().first);  //与最后的type都是3，新的匹配位置往后推
							else path_fea.pop_back();     //最后的type不是3，pop
						}				
						path_fea.push_back(make_pair(one_pair.first, wait_chosen));
						last_index = wait_chosen;
						last_m = one_pair.first;
						//cout << "fixed: " << one_pair.first << "   " << one_pair.second << endl;
					}
					else   //con_mid[one_pair.first].type == fea_p
					{
						int repet_index = -1;
						FOR(pin, 0, path_fea.size())
						{
							if (waiting_merge[min_index] == path_fea[pin].second)
							{
								repet_index = pin;
							}
						}
						if (path_fea.empty() || repet_index == -1)
						{
							path_fea.push_back(make_pair(one_pair.first, waiting_merge[min_index]));
							last_index = waiting_merge[min_index];
							last_m = one_pair.first;
						}
						//repet_index != -1 & same type
						else if (con_mid[one_pair.first].type == con_mid[path_fea[repet_index].first].type)
						{
							int cand_index = (waiting_merge[min_index] + shift) % cand_tar.size();
							int path_rep = path_fea[repet_index].first;
							int wait_min = waiting_merge[min_index];
							double normal_r = 0.5 / (2 * PI);

							double ang1 = angle_2v(con_mid_c[(one_pair.first - 1 + pcsize) % pcsize] - con_mid_c[one_pair.first], con_mid_c[(one_pair.first + 1) % pcsize] - con_mid_c[one_pair.first]);
							double ang2 = angle_2v(con_mid_c[(path_rep - 1 + pcsize) % pcsize] - con_mid_c[path_rep], con_mid_c[(path_rep + 1) % pcsize] - con_mid_c[path_rep]);
							double ang_cand = angle_2v(con2_c[(wait_min - 1 + pcsize) % pcsize] - con2_c[wait_min], con2_c[(wait_min + 1) % pcsize] - con2_c[wait_min]);

							double one_tar = tar_length_2p(inner_tar[one_pair.first], cand_tar[cand_index]) + normal_r*abs(ang1 - ang_cand);
							double back_tar = tar_length_2p(inner_tar[path_rep], cand_tar[cand_index]) + normal_r*abs(ang2 - ang_cand);

							if (one_tar < back_tar)
							{
								path_fea.pop_back();
								path_fea.push_back(make_pair(one_pair.first, waiting_merge[min_index]));
								last_index = waiting_merge[min_index];
								last_m = one_pair.first;
							}
						}
					}
				}
			}
		}
		FOR(mm, 0, path_fea.size()) cout << path_fea[mm].first << "   " << path_fea[mm].second << endl;
		cand_fea_paths.push_back(path_fea);
	}

	vector<pair<int, bool>> Tiling_opt::compare_TAR(vector<Point_f> contour_mid, int chosen_num, int window_width) //chosen_num  选择的最终结果的数目
	{
		//int all_types = 2;
		//std::cout << "inner_c num:" << inner_c.size() << endl;
		vector<Point2f> contour_mid_ = conf_trans(contour_mid);
		double iso_mid = isoperimetric_inequality(contour_mid_);
		cout << "Compare with " << all_types << " types candadites." << endl;

		vector<pair<int, bool>> all_final;
		vector<double> all_result;
		double shape_com_mid;
		vector<vector<double>> tar_mid = compute_TAR(contour_mid_, shape_com_mid);
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
			vector<Point_f> sec = contour_dataset[index];
			vector<vector<double>> tar_sec = all_con_tars[index];
			vector<vector<double>> tar_sec_f = all_con_tars_flip[index];
			vector<pair<int, int>> path;
			int shift = 0;
			double dis_sc = abs(shape_com_mid - all_shape_complexity[index]);
			int cmsize = contour_mid.size();
			vector<int> mid_index(cmsize,0), sec_f(cmsize,0), sec_fflip(cmsize, 0);
			FOR(n, 0, cmsize)
			{
				mid_index[n]= contour_mid[n].type;
				sec_f[n] = sec[n].type;
				sec_fflip[n] = sec[cmsize - 1 - n].type;
				//cout << n << "   " << mid_index[n] << "    " << sec_f[n] << "     " << sec_fflip[n] << endl;
			}
			double re = tar_mismatch_fea(tar_mid, tar_sec, mid_index, sec_f, path, shift, window_width);
			//double re = tar_mismatch(tar_mid, tar_sec, path, shift, window_width);
			//re = re / (1 + shape_com_mid + all_shape_complexity[index]);  //ori
			re = re / (1 + shape_com_mid + all_shape_complexity[index] - 5 * dis_sc - 2 * diss);
			re = 1 - 0.5*re / path.size();
			//cout << "re: " << re << endl;
			double re2 = tar_mismatch_fea(tar_mid, tar_sec_f, mid_index, sec_fflip, path, shift, window_width);
			//double re2 = tar_mismatch(tar_mid, tar_sec_f, path, shift, window_width);
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
		cout << ttrt<< "  patterns beyond ShapeComplexity DisThres: "  << endl;
		return all_total_mid;
	}

	void Tiling_opt::feature_match(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<int> first_fea, vector<int> sec_fea, vector<pair<int, int>>& path, int &sec_shift, int width)
	{
		vector<vector<double>> first_fea_arr;
		vector<vector<double>> second_fea_arr;
		int fsize = first_fea.size();
		int ssize = sec_fea.size();
		FOR(i, 0, fsize)
		{
			first_fea_arr.push_back(first_arr[first_fea[i]]);
		}
		FOR(i, 0, ssize)
		{
			second_fea_arr.push_back(second_arr[sec_fea[i]]);
		}
		double re = tar_mismatch(first_fea_arr, second_fea_arr, path, sec_shift, width);
		vector<int> sec_fea2;
		vector<vector<double>> second_fea_arr2;
		FOR(i, 0, ssize)
		{
			sec_fea2.push_back(sec_fea[(i + sec_shift)% ssize]);
			second_fea_arr2.push_back(second_fea_arr[(i + sec_shift) % ssize]);
		}
		sec_fea = sec_fea2;
		vector<pair<int, int>> waiting_merge;
		vector<pair<int, int>> path_new;
		FOR(i, 0, path.size())
		{
			if(waiting_merge.empty()) 
				waiting_merge.push_back(path[i]);
			else
			{
				if (path[i].first == waiting_merge[0].first || path[i].second == waiting_merge[0].second)
				{
					waiting_merge.push_back(path[i]);
				}
				else
				{
					int min_index = 0;
					double min_dis_tar = 10000;
					FOR(j, 0, waiting_merge.size())
					{
						//calculate the pair with less tar distance
						double dis_tar= tar_length_2p(first_fea_arr[waiting_merge[j].first], second_fea_arr2[waiting_merge[j].second]);
						if (dis_tar < min_dis_tar)
						{
							min_dis_tar = dis_tar;
							min_index = j;
						}
					}
					path_new.push_back(waiting_merge[min_index]);
					waiting_merge.swap(vector<pair<int, int>>());
					waiting_merge.push_back(path[i]);
				}
			}
		}
		//处理waiting_merge中剩下的pair
		int min_index = 0;
		double min_dis_tar = 10000;
		FOR(j, 0, waiting_merge.size())
		{
			//calculate the pair with less tar distance
			double dis_tar = tar_length_2p(first_fea_arr[waiting_merge[j].first], second_fea_arr2[waiting_merge[j].second]);
			if (dis_tar < min_dis_tar)
			{
				min_dis_tar = dis_tar;
				min_index = j;
			}
		}
		path_new.push_back(waiting_merge[min_index]);
		path = path_new;
	}

	vector<Point_f> Tiling_opt::morphing(vector<Point_f> contour1, vector<Point_f> contour2, vector<pair<int, int>> final_pair, double ratio)
	{
		int fpsize = final_pair.size();
		int cnum1 = contour1.size();
		int cnum2 = contour2.size();
		Point2f cen1 = center_p(conf_trans(contour1));
		//cout << "Matched feature point pair : " << fpsize << endl;
		//将两个轮廓分段存储
		vector<vector<Point_f>> contour1_seg;
		vector<vector<Point_f>> contour2_seg;
		vector<Point_f> final_con;
		for (int g = 0; g < fpsize; g++)
		{
			vector<Point_f> c_seg1;
			vector<Point_f> c_seg2;
			//cout << final_pair[g].first << "  :  " << final_pair[g].second << endl;
			int s_index1 = final_pair[g].first;
			int e_index1 = final_pair[(g + 1) % fpsize].first;
			e_index1 > s_index1 ? e_index1 : e_index1 += cnum1;
			//e_index1 = e_index1 + (e_index1 > s_index1 ? 0 : cnum1);
			int s_index2 = final_pair[g].second;
			int e_index2 = final_pair[(g + 1) % fpsize].second;

			e_index2> s_index2 ? e_index2 : e_index2 += cnum2;
			for (int n = s_index1; n <= e_index1; n++)
				c_seg1.push_back(contour1[n%cnum1]);
			for (int n = s_index2; n <= e_index2; n++)
				c_seg2.push_back(contour2[n%cnum2]);
			contour1_seg.push_back(c_seg1);
			contour2_seg.push_back(c_seg2);
		}

		//cout << "contour1_seg:  " << contour1_seg.size() << "    contour2_seg:  " << contour2_seg.size() << endl;
		//	<< contour1_seg.back().back().point << "   " << contour1_seg.back().back().type << endl;
		Mat draw_seg = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift_seg = Point2f(0.5*draw_seg.cols, 0.5*draw_seg.rows) - cen1;
		Point_f start_new = Point_f(ratio*contour1_seg[0][0].point + (1 - ratio)*contour2_seg[0][0].point, contour1_seg[0][0].type);
		double num_e = 0;
		FOR(seg_index, 0, fpsize)
		{
			if (seg_index == fpsize - 1) {
				contour1_seg[seg_index].back().type = 4;
			}
			vector<Point_f> morph_seg = morph_segment(contour1_seg[seg_index], contour2_seg[seg_index], start_new, ratio, num_e);
			final_con.insert(final_con.end(), morph_seg.begin(), morph_seg.end() - 1);
			start_new = morph_seg.back();
			line(draw_seg, morph_seg[0].point + shift_seg, morph_seg.back().point + shift_seg, Scalar(50, 225, 50), 1);
			//cout << contour1_seg[seg_index].size()<<"  "<< contour2_seg[seg_index].size()<<"   "<<morph_seg.size()<<"   "<< morph_seg[0].point << "   " << morph_seg.back().point << endl;
			draw_contour_points(draw_seg, conf_trans(contour1_seg[seg_index]), shift_seg, 5, 2);
			draw_contour_points(draw_seg, conf_trans(contour2_seg[seg_index]), shift_seg+(contour1_seg[seg_index][0].point- contour2_seg[seg_index][0].point), 7, 2);
			
			//cout << morph_seg.size() << "   "<< morph_seg[0].point<< morph_seg[1].point<< morph_seg[2].point<<endl<<final_con.size() << final_con[0].point<< final_con[1].point<<endl;
		}
		//cout << "contour1:  " << cnum1 << "    contour2:  " << cnum2 <<"   final: "<< final_con.size()<< endl;
		//FOR(j, 0, final_con.size()) cout << final_con[j].point << endl;
		string save = SavePath;
		imshow("frame.png", draw_seg);
		//imwrite("frame.png", draw_seg);
		Mat draw_result = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift1 = Point2f(0.5*draw_result.cols, 0.5*draw_result.rows) - cen1;
		int dd = 0;
		FOR(tt, 0, contour1_seg.size())
		{
			FOR(i, 0, contour1_seg[tt].size()) 				
				circle(draw_result, contour1_seg[tt][i].point + shift1, 2, Scalar(255, 0, 0), -1);
			FOR(i, 0, contour2_seg[tt].size())
				circle(draw_result, contour2_seg[tt][i].point + shift1, 2, Scalar(0, 255, 0), -1);
		}
		FOR(i, 0, final_con.size())
			circle(draw_result, final_con[i].point + shift1, 2, Scalar(0, 0, 255), -1);
		imshow("morphing result show", draw_result);

		return final_con;

	}

	vector<Point_f> Tiling_opt::morph_segment(vector<Point_f> seg1, vector<Point_f> seg2, Point_f start, double ratio, double &num_error) //start 是上一个的尾端，是固定的
	{
		Point_f end = Point_f(ratio*seg1.back().point + (1- ratio)*seg2.back().point, seg1.back().type);//seg1.back();  //如果end.type==4, end 不需要变以保证首尾相连
		double length_ave = length_2p(start.point, end.point);
		double angle_ave;
		//cout << "1----start: " << start.point << "     --end: " << end.point <<"  type:"<< end.type<< endl;
		//确定框架的参数 
		if (end.type < 4) // 确定新的end点 
		{
			double angle1 = angle_2v(Point2f(1, 0), seg1.back().point - seg1[0].point);
			double angle_12 = cos_2v(seg1.back().point - seg1[0].point, seg2.back().point - seg2[0].point);
			angle_12 = acos(angle_12);
			if (sin_2v(seg1.back().point - seg1[0].point, seg2.back().point - seg2[0].point) < 0) angle_12 = -angle_12;
			angle_ave = angle1 + (1 - ratio)*angle_12;
			length_ave = ratio*length_2p(seg1.back().point, seg1[0].point) + (1 - ratio)* length_2p(seg2.back().point, seg2[0].point);
			//cout << "angle1: " << angle1/PI*180 << " angle_12: " << angle_12 / PI * 180 << " length: " << length_2p(seg1.back().point, seg1[0].point) << "   " << length_2p(seg2.back().point, seg2[0].point) << "  " << length_ave << endl;
			Point2f vec_fin = Point2f(cos(angle_ave), sin(angle_ave));
			end.point = start.point + length_ave*vec_fin;
			//cout << start.point << "   " << length_ave << "    " << vec_fin << endl;

			//double angle1 = cos_2v(Point2f(1, 0), seg1.back().point - seg1[0].point);
			//angle1 = (angle1 > 0.999) ? 0.999 : (angle1 < -0.999) ? -0.999 : angle1; //防止出现acos(1)的情况，会返回错误值

		}
		//cout << "2----start: " << start.point << "     --end: " << end.point << endl;
		Point2f mid_des = start.point + 0.5*length_ave * vertical_vec(end.point - start.point);
		vector<Point2f> aff_des;
		aff_des.push_back(start.point);
		aff_des.push_back(end.point);
		aff_des.push_back(mid_des);

		int number_new = 0.5*(seg1.size() + seg2.size());
		double diff = 0.5*(seg1.size() + seg2.size()) - number_new;
		if ((diff + num_error)==1.0)
		{
			number_new += 1;
			num_error = 0;
		}
		else if (diff > num_error) num_error = diff;
		//按照number_new进行重采样,这里我们在原样的基础上进行变形
		vector<Point2f> seg_sam1 = sampling_seg(conf_trans(seg1), number_new);
		vector<Point2f> seg_sam2 = sampling_seg(conf_trans(seg2), number_new);
		if (seg_sam1.size() != seg_sam2.size() || seg_sam1.size() != number_new)
		{
			cout << "seg_sam1.size() != seg_sam2.size() != number_new" << endl;
		}
		//construct Affine matrix
		vector<Point2f> aff_seg1;
		aff_seg1.push_back(seg_sam1[0]);
		aff_seg1.push_back(seg_sam1.back());
		aff_seg1.push_back(seg_sam1[0] + 0.5*length_2p(seg_sam1.back(), seg_sam1[0]) * vertical_vec(seg_sam1.back() - seg_sam1[0]));
		Mat M_seg1 = getAffineTransform(aff_seg1, aff_des);
		//cout << "M_seg1: " << M_seg1 << endl;
		cv::transform(seg_sam1, seg_sam1, M_seg1);
		vector<Point2f> aff_seg2;
		aff_seg2.push_back(seg_sam2[0]);
		aff_seg2.push_back(seg_sam2.back());
		aff_seg2.push_back(seg_sam2[0] + 0.5*length_2p(seg_sam2.back(), seg_sam2[0]) * vertical_vec(seg_sam2.back() - seg_sam2[0]));
		Mat M_seg2 = getAffineTransform(aff_seg2, aff_des);
		cv::transform(seg_sam2, seg_sam2, M_seg2);

		vector<Point2f> contour_morphed;
		for (int j = 0; j < number_new; j++)
		{
			contour_morphed.push_back(ratio*seg_sam1[j] + (1 - ratio)*seg_sam2[j]);
		}
		//reassign type of point_f
		vector<Point_f> contour_final_;
		FOR(m, 0, contour_morphed.size())  contour_final_.push_back(p2fea(contour_morphed[m], 0));
		contour_final_[0].type = seg1[0].type;
		contour_final_.back().type = seg1.back().type;

		//cout << "Points num of each seg: " << seg1.size() << "  " << seg2.size() << "   " << contour_final_.size() << endl;
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

	void Tiling_opt::merge_contours(vector<Point_f>& c1, vector<Point_f> &c2, vector<int> &anc1, vector<int> &anc2, double ratio)
	{
		int cnum1 = c1.size();
		if (c2.size() != cnum1)
			cout << "Cannot merge for different size" << endl;
		//将两个轮廓分段存储
		vector<vector<Point_f>> c1_seg;
		vector<vector<Point_f>> c2_seg;
		vector<Point_f> final_con;
		for (int g = 0; g < anc1.size(); g++)
		{
			vector<Point_f> c_seg1;
			int s_index1 = anc1[g];
			int e_index1 = anc1[(g + 1) % anc1.size()];
			e_index1 > s_index1 ? e_index1 : e_index1 += cnum1;
			//e_index1 = e_index1 + (e_index1 > s_index1 ? 0 : cnum1);
			for (int n = s_index1; n <= e_index1; n++)
				c_seg1.push_back(c1[n%cnum1]);
			c1_seg.push_back(c_seg1);
			//cout << g << "  " << c_seg1.size() << "  " << c_seg1[0].point << "  " << c_seg1.back().point << endl;
		}
		for (int g = 0; g < anc2.size(); g++)
		{
			vector<Point_f> c_seg2;
			int s_index2 = anc2[g];
			int e_index2 = anc2[(g + 1) % anc2.size()];
			e_index2> s_index2 ? e_index2 : e_index2 += cnum1;
			for (int n = e_index2; n >= s_index2; n--)
				c_seg2.push_back(c2[n%cnum1]);
			c2_seg.push_back(c_seg2);
			//cout << g << "  " << c_seg2.size() << "  " << c_seg2[0].point << "  " << c_seg2.back().point << endl;
		}
		vector<Point_f> c1_;
		vector<Point_f> c2_;
		vector<vector<Point_f>> c2_segs(4, vector<Point_f>());
		double degree_opt = 0;
		double num_e = 0;
		string file_name = "opt_compare.svg";
		std::ofstream file(file_name);
		file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";

		FOR(i, 0, 4)
		{
			vector<Point_f> each_seg; 
			vector<Point2f> each_seg_;
			int seg_size = c1_seg[i].size();
			//cout << i << "  " << seg_size << "  " << c2_seg[(i + 2) % 4].size() << endl;
			Point2f start1 = c1_seg[i][0].point;
			Point2f start2 = c2_seg[(i + 2) % 4][0].point;
			Point2f sh_ = start1 - start2;
			each_seg_ = merge_segment(conf_trans(c1_seg[i]), conf_trans(c2_seg[(i + 2) % 4]), ratio, num_e);


			FOR(p_index, 0, each_seg_.size() - 1) 
				file << "<line x1=\"" << each_seg_[p_index].x << "\" y1=\"" << each_seg_[p_index].y << "\" x2=\"" << each_seg_[p_index + 1].x << "\" y2=\"" << each_seg_[p_index + 1].y << "\" stroke=\"rgb(120,120,120)\" />\n";
			FOR(p_index, 0, each_seg_.size())
				file << "<circle cx=\"" << each_seg_[p_index].x << "\" cy=\"" << each_seg_[p_index].y << "\" r=\"2\" fill=\"rgb(100,100,100)\" />\n";

		    degree_opt = edge_nd_opt(each_seg_, Clock_order);

			FOR(p_index, 0, each_seg_.size() - 1)
				file << "<line x1=\"" << each_seg_[p_index].x+500 << "\" y1=\"" << each_seg_[p_index].y << "\" x2=\"" << each_seg_[p_index + 1].x + 500 << "\" y2=\"" << each_seg_[p_index + 1].y << "\" stroke=\"rgb(149, 202, 135)\" />\n";
			FOR(p_index, 0, each_seg_.size())
				file << "<circle cx=\"" << each_seg_[p_index].x + 500 << "\" cy=\"" << each_seg_[p_index].y << "\" r=\"2\" fill=\"rgb(250, 65, 65)\" />\n";

			cout << "After degree_opt: " << degree_opt << endl;
			for (auto p : each_seg_) each_seg.push_back(Point_f(p, general_p));
			each_seg[0].type = fixed_p;
			each_seg.back().type = fixed_p;

			FOR(m, 0, each_seg.size() - 1)
			{
				c1_.push_back(each_seg[m]);
				Point_f c2_p = each_seg[each_seg.size() - 1 - m];
				c2_p.point = c2_p.point - sh_;
				c2_segs[(i + 2) % 4].push_back(c2_p);
				/*Point2f c2_p = each_seg[each_seg.size() - 1 - m].point - sh_;
				c2_segs[(i + 2) % 4].push_back(Point_f(c2_p, general_p));*/
			}
		}
		file << "</svg>";
		file.close();
		FOR(n, 0, 4)
		{
			for (auto p : c2_segs[n])
				c2_.push_back(p);
		}
		vector<int> anc1_, anc2_;
		FOR(n, 0, c1.size())
		{
			if (c1_[n].type == fixed_p)
				anc1_.push_back(n);
			if (c2_[n].type == fixed_p)
				anc2_.push_back(n);
		}
		anc1 = anc1_;
		anc2 = anc2_;
		c1 = c1_;
		c2 = c2_;
		//cout << c1.size() << "  " << c2.size() << endl;
		Mat merge_ima = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
		Point2f tttsh = Point2f(600, 600) - center_p(conf_trans(c1));
		draw_contour_points(merge_ima, conf_trans(c1), tttsh);
		draw_contour_points(merge_ima, conf_trans(c2), tttsh, 4);
		FOR(cc, 0, 4)
		{
			circle(merge_ima, c1[anc1[cc]].point + tttsh, 2, Scalar(0, 255, 0));
			circle(merge_ima, c2[anc2[cc]].point + tttsh, 2, Scalar(0, 255, 0));
			//cout << c1[anc1[cc]].point << "  c2: " << c2[anc2[cc]].point << endl;
		}
		imshow("After merging: ", merge_ima);
	}

	vector<Point2f> Tiling_opt::merge_segment(vector<Point2f> seg1, vector<Point2f> seg2, double ratio, double &num_error) //以seg1为基准
	{
		Point2f start = seg1[0];
		Point2f end = seg1.back();
		//align seg1 and seg2
		Point2f shfit = start - seg2[0];
		FOR(i, 0, seg2.size()) seg2[i] += shfit;
		int number_new = 0.5*(seg1.size() + seg2.size());
		double diff = 0.5*(seg1.size() + seg2.size()) - number_new;
		if ((diff + num_error) == 1.0)
		{
			number_new += 1;
			num_error = 0;
		}
		else if (diff > num_error) num_error = diff;
		//按照number_new进行重采样,这里我们在原样的基础上进行变形
		vector<Point2f> seg_sam1 = sampling_seg(seg1, number_new);
		vector<Point2f> seg_sam2 = sampling_seg(seg2, number_new);
		if (seg_sam1.size() != seg_sam2.size() || seg_sam1.size() != number_new)
		{
			cout << "seg_sam1.size() != seg_sam2.size() != number_new" << endl;
		}
		vector<Point2f> contour_merged;
		for (int j = 0; j < number_new; j++)
		{
			contour_merged.push_back(ratio*seg_sam1[j] + (1 - ratio)*seg_sam2[j]);
		}

		return contour_merged;

	}

	double Tiling_opt::deform_evalue(vector<Point2f> con, vector<Point2f> ori_c)
	{
		//shape complexity
		double sc, sc_ori;
		vector<vector<double>> tar_c = compute_TAR(con, sc);
		vector<vector<double>> tar_ori = compute_TAR(ori_c, sc_ori);
		double dis_sc = abs(sc - sc_ori);
		//isoperimetric inequality: the ratio of circumference to area
		double iso_c = isoperimetric_inequality(con);
		double iso_ori = isoperimetric_inequality(ori_c);
		double dis_iso = abs(iso_c - iso_ori) / (0.5*iso_c +0.5* iso_ori);
		//cout << "dis_iso: " << dis_iso << endl;
		vector<pair<int, int>> path;
		int shift = 0;
		int csize = con.size();
		double re = tar_mismatch(tar_c, tar_ori, path, shift, WindowsWidth);
		//re = re / (1 + shape_com_mid + all_shape_complexity[index]);  //ori
		re = re / (1 + sc + sc_ori - 5 * dis_sc - 2 * dis_iso);
		re = 0.5*re / path.size();
		//re = 1 - 0.5*re / path.size();
		//cout << "tar_score: " << re << "    ";

		//以下计算面积：T并S-T交S，之前两个ori_c已经对齐，所以可以直接计算, the lower the better
		Point2f sft= center_p(con) - center_p(ori_c);
		FOR(i, 0, csize) ori_c[i] += sft;
		double area_score = evaluation_area(con, ori_c);
		//cout << "area_score: " << area_score << "  ";// << evaluation_area_pixels(con, ori_c) << "    ";

		//计算惩罚项
		double penalty_score = 0; //如果contour1有交叉，会加分;contour2肯定不会有交叉
		int inter_i, inter_j;
		if (self_intersect(con, inter_i, inter_j))
		{
			int margin_l = min(abs(inter_j - inter_i), abs(inter_i + csize - inter_i));
			penalty_score = margin_l / csize;
		}

		double total_score = re + area_score + penalty_score;
		//std::cout << "penalty_score: " << penalty_score <<"    ";
		cout<<"total_score: " << total_score << endl;

		return  total_score;
	}

}



//void Tiling_opt::tiliing_gen_specify2(string nameid, vector<int> anc_points)
//{
//	load_dataset(false);
//	string filepath = DefaultPath;
//	filepath = filepath + "contour/" + nameid + ".txt";
//	prototile_first = protoTile(filepath);
//	string savepath = SavePath;
//	savepath = savepath.substr(0, savepath.length() - 7) + "result_spe/";
//	savepath += nameid;
//	const char *na = savepath.c_str();
//	if (_access(na, 0) != -1) printf("The  file/dir had been Exisit \n");
//	else	_mkdir(na);
//	vector<int> cand_points = prototile_first.getCandPoints(prototile_first.contour_f);
//	int trans = Tanslation_rule_spec(cand_points, prototile_first.contour_f, savepath, anc_points);
//	//int flips = Flipping_rule(p_p_index, cont_orig, rootname);
//	int all_inner_num = all_inner_conts.size();
//	if (all_inner_num == 0)
//	{
//		std::cout << "no right placement" << endl;
//		return;
//	}
//	FOR(i, 0, 1)
//		//FOR(i, 0, all_inner_num)
//	{
//		match_candidate_(i);
//		//feature_match();
//		for (int j = 0; j < 1; j++)
//		{
//			vector<pair<int, int>> path = cand_paths[j];
//			cout << "path num: " << path.size() << endl;
//			vector<Point_f> contour1 = prototile_mid.contour_f;
//			vector<Point_f> contour2 = candidate_contours[j];
//			prototile_second = protoTile(contour2);
//			vector<pair<int, int>> path_fea = cand_fea_paths[j];
//
//			//show the path
//			Point2f sh = Point2f(300, 300) - center_p(conf_trans(contour1));
//			Point2f shh2 = sh + Point2f(400, 0);
//			Mat fea_match = Mat(600, 1000, CV_8UC3, Scalar(255, 255, 255));
//			draw_contour_points(fea_match, prototile_mid.contour, sh, 5, 2);
//			draw_contour_points(fea_match, prototile_second.contour, shh2, 7, 2);
//			FOR(i, 0, prototile_mid.feature_points.size())
//				circle(fea_match, prototile_mid.contour[prototile_mid.feature_points[i]] + sh, 3, Scalar(0, 0, 255), -1);
//			FOR(i, 0, prototile_second.feature_points.size())
//				circle(fea_match, prototile_second.contour[prototile_second.feature_points[i]] + shh2, 3, Scalar(0, 125, 255), -1);
//			circle(fea_match, prototile_mid.contour[prototile_mid.feature_points[0]] + sh, 5, Scalar(120, 0, 255), -1);
//			circle(fea_match, prototile_second.contour[prototile_second.feature_points[0]] + shh2, 5, Scalar(0, 125, 255), -1);
//			circle(fea_match, prototile_mid.contour[0] + sh, 6, Scalar(25, 200, 25), 2);
//			circle(fea_match, prototile_second.contour[0] + shh2, 6, Scalar(25, 200, 25), 2);
//			FOR(i, 0, path_fea.size())
//			{
//				//int tsfsize = prototile_second.feature_points.size();
//				Point2f f1 = prototile_mid.contour[path_fea[i].first] + sh;
//				Point2f f2 = prototile_second.contour[path_fea[i].second] + shh2;
//				line(fea_match, f1, f2, colorbar[6].second);
//			}
//			imshow("fea match", fea_match);
//			//----------show the feature------------
//
//			double con_sc;
//			vector<vector<double>> contour2_tar;
//			//求contour_对应的point_f和tar值
//			//contour2_tar = computeTAR(contour_, con_sc, 0.5);
//			double angle = 165;
//			double ratio_max = 100;
//			double result_score = 0;
//			double ratio = 0.5;
//			vector<Point_f> contour_2 = morphing(contour1, contour2, path_fea, ratio);
//			cout << "contour_2.size:  " << contour_2.size() << "   " << contour1.size() << endl;
//			//FOR(mm, 0, contour_2.size()) cout << contour_2[mm].type << "   " << contour1[mm].type << endl;
//			//vector<Point_f> contour_2 = morphing_dir(contour1, contour2, path, ratio);
//			Mat pair_match = Mat(600, 1000, CV_8UC3, Scalar(255, 255, 255));
//			draw_pair(pair_match, conf_trans(contour1), conf_trans(contour2), path, sh);
//			Point2f sh2 = Point2f(700, 300) - center_p(conf_trans(contour_2));
//			draw_contour_points(pair_match, conf_trans(contour_2), sh2, 3, 2);
//			imshow("correspond path", pair_match);
//
//			//calculate two deformed shapes
//			vector<Point_f> con_re;
//			vector<int> anc_mid;
//			vector<int> anc_re;
//			int csize = contour_2.size();
//			FOR(t, 0, csize)
//			{
//				if (contour_2[t].type == fixed_p)
//					anc_mid.push_back(t);
//			}
//			Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
//			if (translation_spec(contour_2, con_re, anc_mid, anc_re, draw2))
//			{
//				cout << "OK! No intersection!" << endl;
//			}
//
//			vector<Point2f> contour_dst = conf_trans(contour_2);
//			vector<Point2f> cont_re = conf_trans(con_re);
//			Mat two_c = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
//			draw_contour_points(two_c, contour_dst, Point2f(400, 400) - center_p(contour_dst), 5, 2);
//			draw_contour_points(two_c, cont_re, Point2f(400, 400) - center_p(contour_dst), 6, 2);
//			imshow("2 contours:", two_c);
//			MatrixXd V;
//			MatrixXi F;
//			//vector<Point2f> tt = triangulate_2Contours(contour_dst, cont_re,V,F);
//			//vector<Point2f> tt = triangulate_Contours_bbx(cont_re, anc_re, V, F);
//			//vector<Point2f> tt = triangulate_Contours_bbx(contour_dst, anc_mid, V, F);
//			vector<Point2f> tt = triangulate_bbx(contour_dst, V, F);
//			Mat image = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
//			Point stf(0, 600);
//			for (size_t i = 0; i < F.rows(); i++)
//			{
//				Point pt1(V.row(F(i, 0)).x(), V.row(F(i, 0)).y());
//				Point pt2(V.row(F(i, 1)).x(), V.row(F(i, 1)).y());
//				Point pt3(V.row(F(i, 2)).x(), V.row(F(i, 2)).y());
//				//cout << F(i, 0) << "   " << F(i, 1) << "   " << F(i, 2) << endl;
//				//cout << pt1 << "    " << pt2 << "   " << pt3 << endl;
//				line(image, pt1 + stf, pt2 + stf, Scalar(0, 255, 0), 1, LINE_AA);
//				line(image, pt2 + stf, pt3 + stf, Scalar(0, 255, 0), 1, LINE_AA);
//				line(image, pt3 + stf, pt1 + stf, Scalar(0, 255, 0), 1, LINE_AA);
//			}
//			//imwrite("Triangulation.png", image);
//			imshow("Triangulation", image);
//			double degree_after_opt = 0;
//			if (pers_trans)
//			{
//				vector<Point2f> frame = { contour_dst[anc_mid[0]], contour_dst[anc_mid[1]], contour_dst[anc_mid[2]], contour_dst[anc_mid[3]] };
//				//vector<Point2f> frame = { cont_re[anc_re[0]], cont_re[anc_re[1]], cont_re[anc_re[2]], cont_re[anc_re[3]] };
//				int min_type = 0;
//				double min_l = 10000;
//				vector<Point2f> frame_b;// = base_frame(frame, 0);
//				FOR(f_type, 0, 4)
//				{
//					vector<Point2f> frame_bb = base_frame(frame, f_type);
//					vector<pair<int, int>> path_min = { make_pair(0,0),make_pair(1,1),make_pair(2,2),make_pair(3,3) };
//					double alige_e = conotour_align(frame, frame_bb, path_min);
//					if (alige_e < min_l)
//					{
//						min_l = alige_e;
//						min_type = f_type;
//						frame_b = frame_bb;
//					}
//				}
//				cout << "frame_type: " << frame_type[min_type] << "   " << min_l << endl;
//				vector<int> anc_mid_;
//				FOR(g, 0, 4)
//				{
//					cout << "frame[g]: " << frame[g] << endl;
//					anc_mid_.push_back(point_locate(tt, frame[g]));
//				}
//
//				// 构建handle area
//				//vector<int> handle_area;
//				//vector<Point2f> handle_points;
//				//FOR(index_num, 0, anc_mid.size())
//				//{
//				//	Point2f shift_new = frame_b[index_num] - frame[index_num];
//				//	int area_width = handle_area_width;
//				//	vector<Point2f> handle_one;
//				//	FOR(width_n, -area_width, area_width + 1) 
//				//	//for (int width_n = -area_width; width_n < area_width + 1; width_n += 2)
//				//	{
//				//		int index_ = (anc_mid[index_num] + width_n + csize) % csize;
//				//		handle_area.push_back(index_);
//				//		//handle_points.push_back(contour_2[index_].point + shift_new);
//				//		handle_one.push_back(contour_2[index_].point + shift_new);
//				//	}
//				//	//cout << "handle_one before: " << handle_one[0] << "  " << handle_one[1] << "  " << handle_one[1] << endl;
//				//	int cen_index = handle_one.size() / 2;
//				//	Point2f cen_line = handle_one[0] + handle_one[handle_one.size() - 1] - 2*handle_one[cen_index];
//				//	double cos1 = cos_2v(cen_line, shift_new);
//				//	double sin1 = sin_2v(cen_line, shift_new);
//				//	double angle = acos(cos1) / PI * 180;
//				//	if (sin1 < 0) angle = -angle;		
//				//	if (abs(angle)>90) angle = 0.5*(angle - angle / abs(angle) * 180);
//				//	else angle = 0.5*angle;
//				//	//cout << cen_line << " " << shift_new << cos1 << "  " << angle << endl;
//				//	handle_one = Rotate_contour(handle_one, handle_one[cen_index], angle);
//				//	//cout << "handle_one after: " << handle_one[0] << "  " << handle_one[1] << "  " << handle_one[1] << endl;
//				//	for(auto p: handle_one)
//				//		handle_points.push_back(p);
//				//}
//				Mat draw_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
//				Point2f sh1 = Point2f(400, 400) - center_p(cont_re);
//				draw_contour(draw_, cont_re, sh1, 5);
//				FOR(m, 0, 4) circle(draw_, frame[m] + sh1, 3, Scalar(125, 0, 0));
//				//
//				//MatrixXd V;
//				//MatrixXi F;
//				//vector<Point2f> con = contour_dst;
//				////fileout("D:/vs2015project/Dihedral_new/Dihedral_new/mid_shape.txt", con);
//				//int consize = con.size();
//				//triangulateContour(con, V, F);
//
//				//vector<Point2f> ttt = { tt[0],tt[1],tt[2],tt[3] };
//				//Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
//				//perspectiveTransform(ttt, frame_b, rot_mat);
//				/*vector<int> amid_ = { anc_mid_[0] + 4,anc_mid_[1] + 4,anc_mid_[2] + 4,anc_mid_[3] + 4 };
//				anc_mid_ = amid_;*/
//
//				string parap = ParaPath;
//				string obj_path = parap + "mid_result.obj";
//				string para_path = parap + "deform_para.txt";
//				string deformed_c = parap + "deformed_c.txt";
//				write_obj(obj_path, V, F);
//				write_para(para_path, anc_mid_, frame_b);
//				//write_para(para_path, anc_mid, frame_b);
//				//write_para(para_path, handle_area, handle_points);
//				string command = "D:/vs2015project/Dihedral_new/Dihedral_new/ARAP_Deform.exe  " + obj_path + "  " + para_path + " " + deformed_c;
//				//string command = "D:/vs2015project/ARAP_Deform/x64/Debug/ARAP_Deform.exe  " + obj_path + "  " + para_path + " " + deformed_c;
//				cout << command << endl;
//				system(command.c_str());
//				contour_dst = load_point_file(deformed_c);
//				vector<Point2f> con_tem;
//				FOR(m, 0, csize) con_tem.push_back(contour_dst[m]);
//				contour_dst = con_tem;
//
//				draw_contour(draw_, contour_dst, sh1, 8);
//				FOR(m, 0, 4) circle(draw_, frame_b[m] + sh1, 3, Scalar(0, 255, 0));
//				imshow("After align handle points:", draw_);
//				con_re = set_flags(contour_dst, con_re);
//				//Mat draw_ = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
//				//vector<Point2f> frame = { contour_2[anc_mid[0]].point, contour_2[anc_mid[1]].point, contour_2[anc_mid[2]].point, contour_2[anc_mid[3]].point };
//				//vector<Point2f> frame_b = base_frame(frame, 0);
//				//Point2f sh1(0, 400);
//				//draw_contour(draw_, contour_dst, sh1, 5);
//				//FOR(m, 0, 4) circle(draw_, frame[m] + sh1, 3, Scalar(125, 0, 0));
//				//Mat rot_mat = getPerspectiveTransform(frame, frame_b);// getAffineTransform(frame_1, frame_b_1);
//				//perspectiveTransform(contour_dst, contour_dst, rot_mat);
//				//draw_contour(draw_, contour_dst, sh1, 8);
//				//FOR(m, 0, 4) circle(draw_, frame_b[m] + sh1, 3, Scalar(0, 255, 0));
//				//imshow("123123", draw_);
//			}
//			/*contour_2 = set_flags(contour_dst, contour_2);
//			Mat draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
//			int first = 0, second = 0;
//			if (self_intersect(contour_dst, first, second)) cout << "Bad  intersection!" << endl;
//			if (translation_spec(contour_2, con_re, anc_mid, anc_re, draw2))
//			{
//			cout << "OK! No intersection!" << endl;
//			}
//			contour_dst = conf_trans(con_re);*/
//			//con_re = set_flags(contour_dst, con_re);
//			if (coll_opt)
//			{
//				//contour_de_crossing(contour_dst);
//				contour_fine_tuning(contour_dst);
//				con_re = set_flags(contour_dst, con_re);
//			}
//			if (deve_opt)
//			{
//				vector<Point2f> resam_dst = sampling_ave(contour_dst, contour_dst.size());
//				vector<Point2f> resam_ = { contour_dst[anc_re[0]],contour_dst[anc_re[1]], contour_dst[anc_re[2]], contour_dst[anc_re[3]] };
//				anc_re = relocate(resam_, resam_dst);
//				FOR(ii, 0, 4)
//				{
//					resam_dst[anc_re[ii]] = resam_[ii];
//					cout << resam_[ii] << "   " << resam_dst[anc_re[ii]] << endl;
//				}
//
//				degree_after_opt = whole_con_opt(resam_dst, anc_re, 0);
//				//degree_after_opt = whole_con_opt(contour_dst, anc_re, 0);
//				//degree_after_opt = whole_con_opt(contour_dst, anc_mid, 0);
//				cout << "After developable optimation, the collision degree: " << degree_after_opt << endl;
//				vector<Point_f> con_resam;
//				FOR(ii, 0, resam_dst.size()) con_resam.push_back(Point_f(resam_dst[ii], general_p));
//				FOR(jj, 0, anc_re.size()) con_resam[anc_re[jj]].type = fixed_p;
//				con_re = con_resam;
//			}
//
//			draw2 = Mat(1200, 1600, CV_8UC3, Scalar(255, 255, 255));
//			if (translation_spec(con_re, contour_2, anc_re, anc_mid, draw2))
//			{
//				cout << "OK! No intersection!" << endl;
//			}
//			else cout << "Bad result with intersection!" << endl;
//			protoTile c1, c2;
//			c1.set_contour(conf_trans(con_re), anc_re);
//			c2.set_contour(conf_trans(contour_2), anc_mid);
//			RotationVis(c1, c2, AntiClockWise);
//
//		}
//	}
//}